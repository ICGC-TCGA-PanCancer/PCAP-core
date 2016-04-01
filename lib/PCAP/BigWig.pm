package PCAP::BigWig;

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014-2016 ICGC PanCancer Project
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

use strict;
use autodie qw(:all);
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use File::Spec;

use PCAP::Threaded;

sub bamToBw {
  my ($index, $options) = @_;

  my $tmp = $options->{'tmp'};

  return if(exists $options->{'index'} && $index != $options->{'index'});
  return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @seqs = @{$options->{'sequences'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    my $outfile = File::Spec->catfile($options->{'tmp'}, $seq.'.bw');

    my $command = q{bash -c 'set pipefail; };
    if($options->{'bam'} =~ m/\.bam$/) {
      $command .= _which('bam2bedgraph');
      $command .= q{ }.$options->{'bam'}.q{ }.$seq;
    }
    else {
      $command .= _which('samtools');
      $command .= q{ view -T }.$options->{'reference'};
      $command .= q{ -ub }.$options->{'bam'}.q{ }.$seq;
      $command .= ' | '._which('bam2bedgraph').' - ';
    }
    $command .= ' | ';
    $command .= _which('wigToBigWig');
    $command .= ' -fixedSummaries -keepAllChromosomes stdin '.$options->{'reference'}.'.fai '.$outfile;
    $command .= q{'};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub mergeBw {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my @files;
  opendir(my $dh, $options->{'tmp'});
  while(readdir $dh) {
    next if($_ =~ m/^[.]/);
    push @files, File::Spec->catfile($options->{'tmp'}, $_) if($_ =~ m/[.]bw$/);
  }

  my $bedGraph = File::Spec->catfile($options->{'tmp'}, 'merged.bedGraph');
  my $command = _which('bigWigMerge').q{ };
  $command .= join q{ }, @files;
  $command .= q{ }.$bedGraph;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub generateBw {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $bedGraph = File::Spec->catfile($options->{'tmp'}, 'merged.bedGraph');

  my $outfile = File::Spec->catfile($options->{'outdir'}, 'merged.bw');
  my $command = _which('wigToBigWig');
  $command .= ' '.$bedGraph.' '.$options->{'reference'}.'.fai '.$outfile;
  #wigToBigWig /var/tmp/kr2/bwTest/X.bedGraph /var/tmp/kr2/bwTest/X.bedGraph /var/tmp/kr2/bwTest/new.bw

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

1;

__END__

=head1 PCAP::BigWig

Support module to generate BigWig coverage files from a BAM or CRAM file.

=head2 METHODS

=over 2

=item bamToBw

Generates BigWig files on a pre-chromosome basis to allow parallel (or short-recovery) processing.

=item mergeBw

Merges the full set of BigWig files into an intermediate BedGraph file

=item generateBw

Converts merged BedGraph to a full BigWig file.

=back
