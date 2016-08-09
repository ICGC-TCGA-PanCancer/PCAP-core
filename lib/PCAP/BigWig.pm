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

  my $filter = 3844; # see https://broadinstitute.github.io/picard/explain-flags.html
  $filter = $options->{'filter'} if(exists($options->{'filter'});

  my @seqs = @{$options->{'sequences'}};
  my $iter = 1;
  for my $seq(@seqs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    my $outfile = q{'}.File::Spec->catfile($options->{'tmp'}, $seq.'.bw').q{'};

    my $command = q{bash -c "set pipefail; };
    $command .= _which('bam2bedgraph');
    $command .= q{ -f }.$filter;
    $command .= q{ -r }.$seq;
    $command .= q{ -i }.$options->{'bam'};
    $command .= ' | ';
    $command .= _which('wigToBigWig');
    $command .= ' -fixedSummaries -keepAllChromosomes stdin '.$options->{'reference'}.'.fai '.$outfile;
    $command .= q{"};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub generateBw {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $outfile = File::Spec->catfile($options->{'outdir'}, $options->{'sample'}.'.bw');

  my $command = sprintf '%s -p %s -f %s -o %s', _which('bwcat'),
                                                $options->{'tmp'},
                                                $options->{'reference'}.'.fai',
                                                $outfile;

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
