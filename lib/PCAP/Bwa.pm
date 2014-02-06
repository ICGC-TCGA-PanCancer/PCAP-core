package PCAP::Bwa;

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
use File::Spec;
use File::Which qw(which);

use PCAP::Bwa::Meta;

const my $BWA_ALN => q{ aln%s -t %s -f %s_%s.sai %s %s.%s};
const my $BAMFASTQ => q{ exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=%s filename=%s};
const my $BWA_MEM => q{ mem%s -M -T 0 -R %s -t %s %s};
const my $ALN_TO_SORTED => q{ sampe -P -a 1000 -r '%s' %s %s_1.sai %s_2.sai %s.%s %s.%s | %s fixmate=1 inputformat=sam level=1 tmpfile=%s_tmp O=%s_sorted.bam};
const my $BAMSORT => q{ fixmate=1 inputformat=sam level=1 tmpfile=%s_tmp O=%s_sorted.bam inputthreads=%s outputthreads=%s};

sub bwa_mem {
  # uncoverable subroutine
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my $input_meta = $options->{'meta_set'};
  my $index = 0;
  for my $input(@{$input_meta}) {
    $index++;
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
    my $this_stub = $input->tstub;

    my $rg_line;
    if($input->fastq) {
      $rg_line = q{'}.$input->rg_header(q{\t}).q{'};
    }
    else {
      $rg_line = PCAP::Bam::rg_line_for_output($input->in);
    }

    my $bwa = which('bwa') || die "Unable to find 'bwa' in path";

    my $command;

    # uncoverable branch true
    # uncoverable branch false
    if($input->fastq) {
      my $interleaved_fq = q{};
      $interleaved_fq = q{ -p}, unless($input->paired_fq);
      $bwa .= sprintf $BWA_MEM, $interleaved_fq, $rg_line, $options->{'threads'}, $options->{'reference'};
      if($input->paired_fq) {
        $bwa .= ' '.$input->in.'_1.'.$input->fastq; # add correct file extension for various types
        $bwa .= ' '.$input->in.'_2.'.$input->fastq;
      }
      else {
        $bwa .= ' '.$input->in.'.'.$input->fastq;
      }
      $command = $bwa;
    }
    else {
      my $bam2fq = which('bamtofastq') || die "Unable to find 'bwa' in path";
      $bam2fq .= sprintf $BAMFASTQ, File::Spec->catfile($tmp, "bamtofastq.$index"), $input->in;
      $bwa .= sprintf $BWA_MEM, ' -p', $rg_line, $options->{'threads'}, $options->{'reference'}, q{ -};
      $command = "$bam2fq | $bwa";
    }

    my $helpers = 1;
    $helpers = 2 if($options->{'threads'} > 3);

    my $sort = which('bamsort') || die "Unable to find 'bamsort' in path\n";
    $sort .= sprintf $BAMSORT, File::Spec->catfile($tmp, "bamsort.$index"), $input->tstub, $helpers, $helpers;

    $command .= " | $sort";

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

# will use thread count to speed up align but each set of fastq are processed in series
sub bwa_aln {
  # uncoverable subroutine
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my $input_meta = $options->{'meta_set'};
  my $input_counter = 0;
  my $index_counter = 0;
  for my $input(@{$input_meta}) {
    $input_counter++;
    for my $end(1..2) {
      $index_counter++;
      next if(defined $options->{'index'} && $index_counter != $options->{'index'});
      # uncoverable branch true
      # uncoverable branch false
      next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $input_counter, $end);
      my $this_stub = $input->tstub;
      my $command = which('bwa') || die "Unable to find 'bwa' in path";

      # uncoverable branch true
      # uncoverable branch false
      if($input->fastq) {
        $command .= sprintf $BWA_ALN, q{}, $options->{'threads'}, $this_stub, $end, $options->{'reference'}, $input->in.'_'.$end, $input->fastq;
      }
      else {
        $command .= sprintf $BWA_ALN, ' -b -'.$end, $options->{'threads'}, $this_stub, $end, $options->{'reference'}, $this_stub, 'bam';
      }
      PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $input_counter, $end);
      PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $input_counter, $end);
    }
  }
  return 1;
}

sub sampe {
  # uncoverable subroutine
  my ($index, $options) = @_;
  my $tmp = $options->{'tmp'};
  my $input_meta = $options->{'meta_set'};
  my $ref = $options->{'reference'};

  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  # uncoverable branch true
  # uncoverable branch false
  return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my $pathstub = $input_meta->[$index-1]->tstub;
  my $fastq = $input_meta->[$index-1]->fastq;

  my $command = which('bwa') || die "Unable to find 'bwa' in path";;
  my $bamsort = which('bamsort') || die "Unable to find 'bwa' in bamsort";;
  # uncoverable branch true
  # uncoverable branch false
  if(defined $fastq) {
    $command .= sprintf $ALN_TO_SORTED, $input_meta->[$index-1]->rg_header(q{\t}),
                                      $ref,
                                      $pathstub,
                                      $pathstub,
                                      $input_meta->[$index-1]->in.'_1', $fastq,
                                      $input_meta->[$index-1]->in.'_2', $fastq,
                                      $bamsort,
                                      $pathstub,
                                      $pathstub;
  }
  else {
    $command .= sprintf $ALN_TO_SORTED, $input_meta->[$index-1]->rg_header(q{\t}),
                                      $ref,
                                      $pathstub,
                                      $pathstub,
                                      $pathstub, 'bam',
                                      $pathstub, 'bam',
                                      $bamsort,
                                      $pathstub,
                                      $pathstub;
  };
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'sampe', $index);
}

1;

__END__

=head1 NAME

PCAP::Bwa - Generate BWA mappings

=head2 Methods

=over 4

=item bwa_aln

  PCAP::Bwa::bwa_aln($options);

Runs a 'bwa aln' process for a single end of sequencing data.

  options - Hashref, requires the following entries:
    -tmp        : working/output directory depending on application
    -meta_set   : array reference of objects as generated by PCAP::Bwa::Meta::files_to_meta
    -reference  : Path to reference fa[.gz]
    -threads    : Total threads available to process

See L<PCAP::Bwa::Meta::files_to_meta|PCAP::Bwa::Meta/files_to_meta>.

=item bwa_mem

  PCAP::Bwa::bwa_mem($options);

  options - Hashref, requires the following entries:
    -tmp        : working/output directory depending on application
    -meta_set   : array reference of objects as generated by PCAP::Bwa::Meta::files_to_meta
    -reference  : Path to reference fa[.gz]
    -threads    : Total threads available to process

=item sampe

  Pan::Cancer::Bwa::sampe($index, $options)

Runs a 'bwa sampe' process on the specified element of the meta_set found in $options.

  options - Hashref, requires the following entries:
    -tmp        : working/output directory depending on application
    -meta_set   : array reference of objects as generated by PCAP::Bwa::Meta::files_to_meta
    -reference  : Path to reference fa[.gz]

=back
