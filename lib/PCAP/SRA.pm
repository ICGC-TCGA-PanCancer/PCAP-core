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
use Fcntl qw( :mode );

use File::ShareDir qw(module_dir);

use Data::Dumper;

use PCAP::Bam;

const my @REQUIRED_HEADER_TAGS => qw(ID CN PL LB PI SM PU DT PM);
const my @VALID_SEQ_TYPES => qw(WGS WXS RNA);
const my %ABBREV_TO_SOURCE => ( 'WGS' => {'source' => 'GENOMIC',
                                          'selection' => 'RANDOM'},
                                'WXS' => {'source' => 'GENOMIC',
                                          'selection' => 'Hybrid Selection'},
                                'RNA' => {'source' => 'RNA',
                                          'selection' => 'RANDOM'},);
const my @REQUIRED_FIELDS => qw(participant_id sample_id aliquot_id submitter_participant_id
                                submitter_sample_id submitter_aliquot_id dcc_project_code
                                total_lanes);
const my %CV_MAPPINGS => ('analyte_code'     => { 'file' => 'cv_tables/TCGA/portionAnalyte.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'sample_type'      => { 'file' => 'cv_tables/TCGA/sampleType.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'disease_abbr'     => { 'file' => 'cv_tables/TCGA/diseaseStudy.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'tss_id'           => { 'file' => 'cv_tables/TCGA/tissueSourceSite.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_donor_gender' => { 'file' => 'cv_tables/ICGC/dcc_donor_gender.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_primary_site' => { 'file' => 'cv_tables/ICGC/dcc_primary_site.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_project_code' => { 'file' => 'cv_tables/ICGC/dcc_project_code.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_project_name' => { 'file' => 'cv_tables/ICGC/dcc_project_name.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_specimen_type' => {'file' => 'cv_tables/ICGC/dcc_specimen_type.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          'dcc_tumour_subtype' =>{'file' => 'cv_tables/ICGC/dcc_tumour_subtype.txt',
                                                  'column' => 0,
                                                  'header' => 1},
                          );

sub generate_sample_SRA {
  my ($grouped, $options) = @_;
  my $cv_lookups = create_cv_lookups();
  my @analysis_ids;
  my $base_path = $options->{'outdir'};
  for my $seq_type(keys %{$grouped}) {
    for my $sample(keys %{$grouped->{$seq_type}}) {
      my $total_bams = 0;
      for my $lib_id(keys %{$grouped->{$seq_type}->{$sample}}) {
        $total_bams += scalar @{$grouped->{$seq_type}->{$sample}->{$lib_id}};
      }
      for my $lib_id(keys %{$grouped->{$seq_type}->{$sample}}) {
        for my $bam_ob(@{$grouped->{$seq_type}->{$sample}->{$lib_id}}) {
          info_file_data($bam_ob);
          $bam_ob->{'info'}->{'total_lanes'} = $total_bams;
          validate_info($cv_lookups, $bam_ob);

          my $submission_uuid = &uuid;
          my $submission_path = "$base_path/".$submission_uuid;
          make_path($submission_path);
          my %exps;
          my %runs;
          $exps{$bam_ob->{'exp'}} = $bam_ob unless(exists $exps{$bam_ob->{'exp'}});
          push @{$runs{$bam_ob->{'run'}}}, $bam_ob;

          my $run_xmls = run($bam_ob->{'CN'}, \%runs);
          my $exp_xml = experiment_sets($options->{'study'}, \%exps);

          my $analysis_xml = analysis_xml($bam_ob, $options->{'study'}, $sample);
          open my $XML, '>', "$submission_path/analysis.xml";
          print $XML $analysis_xml;
          close $XML;
          my $run_xml = run_set($bam_ob->{'CN'}, $run_xmls);
          open $XML, '>', "$submission_path/run.xml";
          print $XML $run_xml;
          close $XML;
          open $XML, '>', "$submission_path/experiment.xml";
          print $XML $exp_xml;
          close $XML;

          my ($cleaned_filename, $directories, $suffix) = fileparse($bam_ob->{'file'}, '.bam');
          $cleaned_filename .= '.bam';
          symlink abs_path($bam_ob->{'file'}), "$submission_path/$cleaned_filename";
          push @analysis_ids, $submission_uuid;
        }
      }
    }
  }
  print "## Executing the following will complete the submission/upload process:\n";
  my $full_path = abs_path($base_path);
  my $sra_sh_script = "$full_path/auto_upload.sh";
  open my $SH, '>', $sra_sh_script;
  print $SH bash_script($options->{'gnos'}, $full_path, \@analysis_ids);
  close $SH;
  chmod S_IRUSR|S_IXUSR, $sra_sh_script;
  my $log = $sra_sh_script;
  $log .= '.log';
  print "$sra_sh_script >& $log &\n";
  print "tail -f $log\n";
}

sub create_cv_lookups {
  my $data_path = shift; # only for use by test cases
  my %cv_lookup;
  # default to location of running code, if not present
  # likely that module has been installed
  # so try installed area
  unless(defined $data_path && -e $data_path) {
    $data_path = dirname(abs_path($0)).'/../share';
    $data_path = module_dir('PCAP::SRA') unless(-e "$data_path/cv_tables");
  }
  for my $cv_field(keys %CV_MAPPINGS) {
    my $cv_file = "$data_path/$CV_MAPPINGS{$cv_field}{file}";
    die "ERROR: Unable to find controlled vocabulary file for $cv_field: $cv_file\n" unless(-e $cv_file);
    die "ERROR: Controlled vocabulary file for $cv_field is empty: $cv_file\n" unless(-s _);
    open my $CV_IN, '<', $cv_file;
    my @cv_set;
    while(my $line = <$CV_IN>) {
      next if($CV_MAPPINGS{$cv_field}{'header'} && $INPUT_LINE_NUMBER == 1);
      chomp $line;
      push @cv_set, (split /\t/, $line)[$CV_MAPPINGS{$cv_field}{'column'}];
    }
    close $CV_IN;
    $cv_lookup{$cv_field} = \@cv_set;
  }
  return \%cv_lookup;
}

sub validate_info {
  my ($cv_lookup, $bam_ob) = @_;
  my %info = %{$bam_ob->{'info'}};
  for my $key(keys %info) {
    next unless(exists $cv_lookup->{$key});
    die "CV term '$key' has invalid value '$info{$key}' in $bam_ob->{file}\n" unless(first { $info{$key} eq $_ } @{$cv_lookup->{$key}} );
  }
  for my $req(@REQUIRED_FIELDS) {
    die "Required comment field '$req' is missing" unless(exists $info{$req});
  }
  return 1;
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
    $bam->check_paired;
    $bam->read_group_info(\@REQUIRED_HEADER_TAGS);
    my %bam_detail;
    for my $tag(@REQUIRED_HEADER_TAGS) {
      $bam_detail{$tag} = $bam->single_rg_value($tag);
      if($tag eq 'SM') {
        my $sm = lc $bam_detail{$tag};
        if($sm =~ /^([a-f0-9]{8})([a-f0-9]{4})([a-f0-9]{4})([a-f0-9]{4})([a-f0-9]{12})$/) {
          $sm = join q{-}, ($1,$2,$3,$4,$5);
        }
        $bam_detail{$tag} = $sm;
      }
    }
    $bam_detail{'file'} = $bam->{'bam'};
    $bam_detail{'md5'} = $bam->{'md5'};
    $bam_detail{'comments'} = $bam->comments;
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
  my $bam_ob = shift;
  return sprintf '<RUN data_block_name="%s" read_group_label="%s" refcenter="%s" refname="%s"/>'
                  , $bam_ob->{'LB'}
                  , $bam_ob->{'ID'}
                  , $bam_ob->{'CN'}
                  , $bam_ob->{'run'};
}

sub get_md5_from_file {
  my $file = shift;
  open my $IN, '<', $file;
  my $md5 = <$IN>;
  close $IN;
  chomp $md5;
  $md5 =~ s/\s+.*//;
  return $md5;
}

sub analysis_xml {
  my ($bam_ob, $study_name, $aliquot_id) = @_;
  my $dt = $bam_ob->{'DT'};
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
      <TARGET refcenter="OICR" refname="%s" sra_object_type="SAMPLE"/>
    </TARGETS>
    <DATA_BLOCK>
      <FILES>
        %s
      </FILES>
    </DATA_BLOCK>%s
  </ANALYSIS>
</ANALYSIS_SET>
ANALYSISXML
  return sprintf $analysis_xml
                , $bam_ob->{'CN'}
                , $bam_ob->{'DT'}
                , $study_name
                , analysis_run_xml($bam_ob)
                , $aliquot_id
                , file_xml($bam_ob)
                , analysis_attributes($bam_ob->{'info'});
}

sub analysis_attributes {
  my $info = shift;
  my $attr_xml = q{};
  my @tags = sort keys %{$info};
  if(scalar @tags > 0) {
    my @attributes;
    for my $tag(@tags) {
      push @attributes, _attribute_xml($tag, $info->{$tag});
    }
    $attr_xml = "\n    <ANALYSIS_ATTRIBUTES>\n".
                (join "\n", @attributes).
                "\n    </ANALYSIS_ATTRIBUTES>";
  }
  return $attr_xml;
}

sub _attribute_xml {
  my ($tag, $value) = @_;
  my $attr_xml = <<ATTRXML;
      <ANALYSIS_ATTRIBUTE>
        <TAG>%s</TAG>
        <VALUE>%s</VALUE>
      </ANALYSIS_ATTRIBUTE>
ATTRXML
  chomp $attr_xml;
  return sprintf $attr_xml, $tag, $value;
}

sub experiment_sets {
  my ($study, $exp_set) = @_;
  my $experiment_xml = <<EXP_XML;
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.experiment.xsd?view=co">
%s
</EXPERIMENT_SET>
EXP_XML

  my @experiments;
  for my $exp(keys %{$exp_set}) {
    push @experiments, experiment($study, $exp_set->{$exp});
  }
  return sprintf $experiment_xml, (join '', @experiments);
}

sub experiment {
  my ($study, $bam_ob) = @_;
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

  return sprintf $exp_xml , $bam_ob->{'CN'}
                          , $bam_ob->{'exp'}
                          , $study
                          , sample_descriptor($bam_ob)
                          , $bam_ob->{'LB'} # definition on NCBI is incorrect
                          , $bam_ob->{'type'} # WGS, WXS, RNA-Seq
                          , $ABBREV_TO_SOURCE{$bam_ob->{'type'}}->{'source'} # GENOMIC
                          , $ABBREV_TO_SOURCE{$bam_ob->{'type'}}->{'selection'} # Random, Hybrid selection
                          , $bam_ob->{'PI'}
                          , $bam_ob->{'PM'};
}

sub sample_descriptor {
  my $bam_ob = shift;
  my $local_sample = q{};
  if(exists $bam_ob->{'info'}->{'submitter_aliquot_id'}) {
    $local_sample = qq{<SUBMITTER_ID namespace="$bam_ob->{CN}">$bam_ob->{info}->{submitter_aliquot_id}</SUBMITTER_ID>\n          };
  }
  my $samp_desc = <<SAMPXML;
      <SAMPLE_DESCRIPTOR refcenter="OICR" refname="%s">
        <IDENTIFIERS>
          %s<UUID>%s</UUID>
        </IDENTIFIERS>
      </SAMPLE_DESCRIPTOR>
SAMPXML
  chomp $samp_desc;
  return sprintf $samp_desc, $bam_ob->{'SM'}, $local_sample, $bam_ob->{'SM'};
}

sub info_file_data {
  my ($bam_ob) = @_;
  my $info_file = $bam_ob->{'file'}.'.info';
  if(-e $info_file) {
    open my $IN, '<', $info_file;
    while (my $line = <$IN>) {
      chomp $line;
      next if($line eq q{});
      die "$info_file has more that one entry for $1" if(exists $bam_ob->{'info'}->{$1});
      $bam_ob->{'info'}->{$1} = $2 if($line =~ m/^([^:]+):(.*)$/);
    }
  }
  # also check the bam header
  for my $comment(@{$bam_ob->{'comments'}}) {
    next unless($comment =~ m/^([^:]+):(.*)/);
    my ($key, $value) = ($1, $2);
    die "$bam_ob->{file} has more that one entry for $1 (or entry is also in *.info file)" if(exists $bam_ob->{'info'}->{$1});
    if($value =~ /^([a-fA-F0-9]{8})([a-fA-F0-9]{4})([a-fA-F0-9]{4})([a-fA-F0-9]{4})([a-fA-F0-9]{12})$/) {
      $value = join q{-}, ($1,$2,$3,$4,$5);
    }
    $bam_ob->{'info'}->{$key} = $value
  }
  return 1;
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
  my $xml = <<EXPXML;
    <EXPERIMENT_REF refcenter="%s" refname="%s"/>
    <DATA_BLOCK>
      <FILES>
        %s
      </FILES>
    </DATA_BLOCK>
EXPXML
  chomp $xml;
  return $xml;
}

sub _run_xml {
  my $xml = <<RUNXML;
  <RUN center_name="%s" alias="%s">
%s
  </RUN>
RUNXML
  chomp $xml;
  return $xml;
}

sub bash_script {
  my ($gnos_server, $path, $uuids) = @_;
  my $uuid_str = join q{" "}, @{$uuids};
my $script = <<'BASHSCRIPT';
#!/bin/bash
set -e
set -u
set -o pipefail

submitexp=" OK ";
queryext="All matching objects are in a downloadable state";

submit_needed () {
  if [ -e "$1" ]; then
    catres=`cat $1`
    if [ "$catres" != $submitexp ]; then
      return 0
    fi
  else
    return 0
  fi
  return 1
}

upload_needed () {
  uploadlog="$1/gtupload.log"
  if [ -e "$uploadlog" ]; then
    # check against cgquery
    set +e
    tmpfile="$(mktemp)"
    thing="cgquery -s %s analysis_id=$1"
    $thing >& $tmpfile
    if cat "$tmpfile" | grep -q "$queryext"; then
      rm -f $tmpfile
      return 1
    else
      rm -f $tmpfile
      return 0
    fi
  else # no log file so upload needed
    return 0
  fi
  return 1
}

process_uuids () {
  name=$1[@]
  uuids=("${!name}")

  for i in "${uuids[@]}"; do
    submitlog="$i/cgsubmit.log"
    if submit_needed $submitlog; then
      set -x
      cgsubmit -s %s -o $submitlog -u $i -c $GNOS_PERM > $submitlog.out
      set +x
    else
      echo RESUME MESSAGE: cgsubmit previously successful for $i
    fi
    if upload_needed $i; then
      set -x
      gtupload -v -c $GNOS_PERM -u $i/manifest.xml >> $i/gtupload.log 2>&1
      set +x
    else
      echo RESUME MESSAGE: gtupload previously successful for $i
    fi
  done
}

# change into working dir
workarea="%s"
echo Working directory: $workarea
cd $workarea
ids=( "%s" )

process_uuids ids

echo SUCCESSFULLY COMPLETED

BASHSCRIPT
  return sprintf $script, $gnos_server, $gnos_server, $path, $uuid_str;
}

1;

__END__

=head2 Methods

=over 4

=item parse_input

Takes a list of files and converts into basic bam detail structure.

=item get_md5_from_file

Pulls pre-calculated MD5 from file co-located with BAM as *.md5

=item experiment_sets

Generates the EXPERIMENT_SET component of experiment.xml

=item experiment

Generates the EXPERIMENT block for experiment.xml

=item sample_descriptor

Generates the SAMPLE_DESCRIPTOR block for experiment.xml

=item run_set

Generates the RUN_SET block for run.xml

=item run

Generates the RUN block for run.xml and linked EXPERIMENT_REF block

=item info_file_data

Adds any data to be added to ANALYSIS_ATTRIBUTES into the 'info' component of the bam_object.

Data is pulled from bam header @CO field or alternatively from a *.info file co-located with the *.bam

=item analysis_xml

Generates the full analysis.xml content.

=item analysis_attributes

Generates ANALYSIS_ATTRIBUTES block of analysis.xml

=item analysis_run_xml

Generates the RUN XML element for analysis.xml

=item file_xml

Generates the FILE XML element which is used in both analysis.xml and run.xml.

=item group_bams

Group bam files by seq type, sample, library.

=item validate_seq_type

Checks sequencing type is of expected format/type.

=item uuid

Creates a lower-case uuid (GNOS preference)

=item create_cv_lookups

Loads the pre-determined look up files for CV terms.

=item validate_info

Checks any fields in the info block that match the name of CV terms against the valid values.

=item generate_sample_SRA

Generates a set of submission files for each of the grouped bams.
During processing also validates any fields known to have controlled vocab.

=item analysis_xml

Generates the analysis.xml file content.

Takes list of values in this order

  bam object (with info prepopulated)
  study_name
  aliquot_id from BAM RG header SM tag

=item bash_script

Takes output path and list of submission UUIDs.

Generates a bash script that can be run to complete GNOS upload with resume capabilities.

=back
