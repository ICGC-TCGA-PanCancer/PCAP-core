package PCAP::SRA;

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014 ICGC PanCancer Project
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not see:
#   http://www.gnu.org/licenses/gpl-2.0.html
##########LICENCE##########


use PCAP;
our $VERSION = PCAP->VERSION;

use strict;
use autodie qw(:all);
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Const::Fast qw(const);
use Carp qw(croak);
use List::Util qw(first);
use File::Path qw(make_path);
use File::Basename;
use Cwd 'abs_path';
use Data::UUID;
use Data::Dumper;

use PCAP::Bam;

const my @REQUIRED_HEADER_TAGS => qw(ID CN PL LB PI SM PU DT);
const my @VALID_SEQ_TYPES => qw(WGS WXS RNA);
const my %ABBREV_TO_SOURCE => ( 'WGS' => {'source' => 'GENOMIC',
                                          'selection' => 'Random'},
                                'WXS' => {'source' => 'GENOMIC',
                                          'selection' => 'Hybrid Selection'},
                                'RNA' => {'source' => 'RNA',
                                          'selection' => 'Random'},);

sub generate_sample_SRA {
  my ($grouped, $options) = @_;
  my $base_path = $options->{'outdir'};
  for my $seq_type(keys %{$grouped}) {
    for my $sample(keys %{$grouped->{$seq_type}}) {
      for my $lib_id(keys %{$grouped->{$seq_type}->{$sample}}) {
        my $submission_path = "$base_path/".&uuid;
        make_path($submission_path);

        my %exps;
        my %runs;
        my @all_bams;
        for my $bam_ob(@{$grouped->{$seq_type}->{$sample}->{$lib_id}}) {
          $exps{$bam_ob->{'exp'}} = $bam_ob unless(exists $exps{$bam_ob->{'exp'}});
          push @{$runs{$bam_ob->{'run'}}}, $bam_ob;
          push @all_bams, $bam_ob;
        }

        my $centre = $all_bams[0]->{'CN'};

        my $run_xmls = run($centre, \%runs);
        my $exp_xml = experiment_sets($centre, $options->{'study'}, $sample, \%exps);

        my $analysis_xml = analysis_xml($centre, $options->{'study'}, $sample, \@all_bams);
        open my $XML, '>', "$submission_path/analysis.xml";
        print $XML $analysis_xml;
        close $XML;
        my $run_xml = run_set($centre, $run_xmls);
        open $XML, '>', "$submission_path/run.xml";
        print $XML $run_xml;
        close $XML;
        open $XML, '>', "$submission_path/experiment.xml";
        print $XML $exp_xml;
        close $XML;

        for (@all_bams) {
          my ($cleaned_filename, $directories, $suffix) = fileparse($_->{'file'}, '.bam');
          $cleaned_filename .= '.bam';
          symlink abs_path($_->{'file'}), "$submission_path/$cleaned_filename";
        }
      }
    }
  }
}

sub uuid {
  my $ug= new Data::UUID;
  return lc $ug->create_str();
}

sub validate_seq_type {
  my $seq_in = shift;
  die "ERROR: '$seq_in' is not a recognised sequencing/library type\n" unless(_check_seq_type($seq_in));
  return 1;
}

sub _check_seq_type {
  my $seq_in = shift;
  return first {$seq_in eq $_} @VALID_SEQ_TYPES;
}

sub group_bams {
  my ($bam_obs, $in_seq_type) = @_;
  my %grouped;
  for my $bam_ob(@{$bam_obs}) {
    my $sm = $bam_ob->{'SM'};
    my $lb = $bam_ob->{'LB'};
    my ($run) = $bam_ob->{'PU'} =~ m/^[[:alpha:]]+:([^_]+)_[^#]+/;
    $bam_ob->{'run'} = sprintf '%s:%s', $bam_ob->{'CN'}, $run;
    my ($lib_id) = $lb =~ m/^[[:alpha:]]+:[[:alpha:]]+:(.*)$/;
    $bam_ob->{'exp'} = sprintf '%s:%s', $bam_ob->{'run'}, $lib_id;
    my $seq_type;
    if(defined $in_seq_type) {
      $seq_type = $in_seq_type;
    }
    else {
      ($seq_type) = $lb =~ m/^([^\:]+)/;
      die "Valid library type is not encoded in readgroup LB tag: $lb\n" unless(_check_seq_type($seq_type));
    }
    $bam_ob->{'type'} = $seq_type;
    push @{$grouped{$seq_type}{$sm}{$lib_id}}, $bam_ob;
  }
  return \%grouped;
}

sub parse_input {
  my $files = shift;
  my @bam_obs;
  for my $file(@{$files}) {
    my $bam = PCAP::Bam->new($file);
    $bam->read_group_info(\@REQUIRED_HEADER_TAGS);
    my %bam_detail;
    for my $tag(@REQUIRED_HEADER_TAGS) {
      $bam_detail{$tag} = $bam->single_rg_value($tag);
    }
    $bam_detail{'file'} = $bam->{'bam'};
    $bam_detail{'md5'} = $bam->{'md5'};
    push @bam_obs, \%bam_detail
  }
  return \@bam_obs;
}

sub file_xml {
  my $bam = shift;
  my $md5 = get_md5_from_file($bam->{'file'}.'.md5');
  my ($cleaned_filename, $directories, $suffix) = fileparse($bam->{'file'}, '.bam');
  $cleaned_filename .= '.bam';
  return sprintf '<FILE checksum="%s" checksum_method="MD5" filename="%s" filetype="bam"/>'
                          , $md5
                          , $cleaned_filename;
}

sub analysis_run_xml {
  my $bam = shift;
  return sprintf '<RUN data_block_name="%s" read_group_label="%s" refcenter="WTSI" refname="%s"/>'
                  , $bam->{'LB'}
                  , $bam->{'ID'}
                  , $bam->{'run'};
}

sub get_md5_from_file {
  my $file = shift;
  open my $IN, '<', $file;
  my $md5 = <$IN>;
  close $IN;
  chomp $md5;
  return $md5;
}

sub analysis_xml {
  my ($centre_name, $study_name, $aliquot_id, $files) = @_;
  # if assembly short_name is to be used need to add parameter
  # otherwise need to delete assembly section
  my @tmp_dt;
  my @file_xml;
  my @run_xml;
  for(@{$files}) {
    push @tmp_dt, $_->{'DT'};
    push @file_xml, file_xml($_);
    push @run_xml, analysis_run_xml($_);
  }
  @tmp_dt = sort @tmp_dt;
  my $dt = $tmp_dt[-1];
  my $analysis_xml = <<ANALYSISXML;
<ANALYSIS_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.analysis.xsd?view=co">
  <ANALYSIS center_name="%s" analysis_date="%s" >
    <TITLE></TITLE>
    <STUDY_REF refcenter="OICR" refname="%s"/>
    <DESCRIPTION>NA</DESCRIPTION>
    <ANALYSIS_TYPE>
      <REFERENCE_ALIGNMENT>
        <ASSEMBLY>
          <STANDARD short_name="unaligned"/>
        </ASSEMBLY>
        <RUN_LABELS>
          %s
        </RUN_LABELS>
        <SEQ_LABELS>
          <SEQUENCE accession="NA" data_block_name="NA" seq_label="NA"/>
        </SEQ_LABELS>
        <PROCESSING>
          <DIRECTIVES>
            <alignment_includes_unaligned_reads>true</alignment_includes_unaligned_reads>
            <alignment_marks_duplicate_reads>false</alignment_marks_duplicate_reads>
            <alignment_includes_failed_reads>false</alignment_includes_failed_reads>
          </DIRECTIVES>
          <PIPELINE>
            <PIPE_SECTION>
              <STEP_INDEX>NA</STEP_INDEX>
              <PREV_STEP_INDEX>NA</PREV_STEP_INDEX>
              <PROGRAM>BIOBAMBAM</PROGRAM>
              <VERSION>0.0.120</VERSION>
              <NOTES>VENDOR FAILED READS REMOVED BEFORE UPLOAD</NOTES>
            </PIPE_SECTION>
          </PIPELINE>
        </PROCESSING>
      </REFERENCE_ALIGNMENT>
    </ANALYSIS_TYPE>
    <TARGETS>
      <TARGET refcenter="TCGA" refname="%s" sra_object_type="SAMPLE"/>
    </TARGETS>
    <DATA_BLOCK>
      <FILES>
        %s
      </FILES>
    </DATA_BLOCK>
  </ANALYSIS>
</ANALYSIS_SET>
ANALYSISXML
  return sprintf $analysis_xml, $centre_name, $dt, $study_name, (join "\n", @run_xml), $aliquot_id, (join "\n", @file_xml);
}

sub experiment_sets {
  my ($centre, $study, $sample, $exp_set) = @_;
  my $experiment_xml = <<EXP_XML;
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.experiment.xsd?view=co">
%s
</EXPERIMENT_SET>
EXP_XML

  my @experiments;
  for my $exp(keys %{$exp_set}) {
    push @experiments, experiment($centre, $study, $sample, $exp_set->{$exp});
  }
  return sprintf $experiment_xml, (join '', @experiments);
}

sub experiment {
  my ($centre, $study, $sample, $lib) = @_;
  my $exp_xml = <<EXPXML;
  <EXPERIMENT center_name="%s" alias="%s">
    <STUDY_REF refcenter="OICR" refname="%s"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>NA</DESIGN_DESCRIPTION>
%s
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>%s</LIBRARY_NAME>
        <LIBRARY_STRATEGY>%s</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>%s</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>%s</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <PAIRED NOMINAL_LENGTH="%s"/>
        </LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <ILLUMINA>
        <INSTRUMENT_MODEL>%s</INSTRUMENT_MODEL>
      </ILLUMINA>
    </PLATFORM>
  </EXPERIMENT>
EXPXML
  chomp $exp_xml;

warn "sequencer model artificial\n";

  return sprintf $exp_xml , $centre
                          , $lib->{'exp'}
                          , $study
                          , sample_descriptor($lib)
                          , $lib->{'LB'} # definition on NCBI is incorrect
                          , $lib->{'type'} # WGS, WXS, RNA-Seq
                          , $ABBREV_TO_SOURCE{$lib->{'type'}}->{'source'} # GENOMIC
                          , $ABBREV_TO_SOURCE{$lib->{'type'}}->{'selection'} # Random, Hybrid selection
                          , $lib->{'PI'}
                          , 'Illumina HiSeq 2000';#$lib->{'PM'};
}

sub sample_descriptor {
  my $bam_ob = shift;
  my $info = info_file_data($bam_ob);
  my $local_sample = q{};
warn "submitter_id skipped\n";
#  if(exists $info->{'INTERNAL_SAMPLE'}) {
#    $local_sample = qq{\n          <SUBMITTER_ID>$info->{INTERNAL_SAMPLE}</SUBMITTER_ID>};
#  }
  my $samp_desc = <<SAMPXML;
      <SAMPLE_DESCRIPTOR refcenter="OICR" refname="%s">
        <IDENTIFIERS>
          <UUID>%s</UUID>%s
        </IDENTIFIERS>
      </SAMPLE_DESCRIPTOR>
SAMPXML
  chomp $samp_desc;
  return sprintf $samp_desc, $bam_ob->{'SM'}, $bam_ob->{'SM'}, $local_sample;
}

sub info_file_data {
  my $bam_ob = shift;
  my $info_file = $bam_ob->{'file'}.'.info';
  my %info;
  if(-e $info_file) {
    open my $IN, '<', $info_file;
    while (my $line = <$IN>) {
      chomp $line;
      next if($line eq q{});
      $info{$1} = $2 if($line =~ m/^([^:]+):(.*)$/);
    }
  }
  return \%info;
}

sub run_set {
  my ($centre_name, $runs) = @_;
  my $run_set_xml = <<RUN_SET_XML;
<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.run.xsd?view=co">
%s
</RUN_SET>
RUN_SET_XML
  return sprintf $run_set_xml, (join "\n", @{$runs});
}

sub run {
  my ($centre, $run_sets) = @_;
  my @run_xmls;
  for my $run(keys %{$run_sets}) {
    # collate experiments within run
    my %exps;
    my @exp_xmls;
    for my $bam(@{$run_sets->{$run}}) {
      push @{$exps{$bam->{'exp'}}}, $bam;
    }
    for my $exp(keys %exps) {
      my @file_xmls;
      for my $bam(@{$exps{$exp}}) {
        push @file_xmls, file_xml($bam);
      }
      push @exp_xmls, sprintf &_exp_xml, $centre, $exp, (join qq{\n}, @file_xmls);
    }
    push @run_xmls, sprintf &_run_xml, $centre, $run, (join qq{\n}, @exp_xmls);
  }
  return \@run_xmls;
}

sub _exp_xml {
  return <<EXPXML;
    <EXPERIMENT_REF refcenter="%s" refname="%s"/>
    <DATA_BLOCK>
      <FILES>
%s
      </FILES>
    </DATA_BLOCK>
EXPXML
}

sub _run_xml {
  return <<RUNXML;
  <RUN center_name="%s" alias="%s">
%s
  </RUN>
RUNXML
}

1;

__END__

=head2 Methods

=over 4

=item analysis_xml

  Takes list of values in this order

  center_name from BAM RG header CN tag
  analysis_date from BAM RG header DT tag
  study_name
  ??short_name="TODO:Reference Genome" ## possibly dropped??
  aliquot_id from BAM RG header SM tag
  BAM md5
  bam filename

=back
