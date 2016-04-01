#! /usr/bin/perl

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

use strict;
use warnings FATAL=>'all';
use autodie;
use Carp;
use Const::Fast;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Pod::Usage qw(pod2usage);

use PCAP::Bam::Coverage;

my $opts = parse_options();

my $coverage = PCAP::Bam::Coverage->new($opts);

my $depth_data = $coverage->build_depth;

if(defined($opts->{'out'})){
  my $OUT;
  open($OUT, '>', $opts->{'out'}) or croak("Error opening output file for write :$!");
    print $OUT $depth_data,"\n";
  close($OUT);
}else{
  print $depth_data,"\n";
}

sub parse_options{
my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'v|version' => \$opts{'v'},
					'f|xam_file=s' => \$opts{'xam'},
					'r|target_file=s' => \$opts{'target'},
					'o|output_file=s' => \$opts{'out'},
					't|type=s' => \$opts{'type'},
					) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});
  if(defined $opts{'v'}) {
    print PCAP::Bam::Coverage->VERSION,"\n";
    exit 0;
  }

  pod2usage(-message => "Option 'f|xam_file' bam|cram file required.", -verbose => 0) if(!defined $opts{'xam'} || ! -e $opts{'xam'});
  pod2usage(-message => "Option 'r|target_file' target gff3|bed file required.", -verbose => 0) if(!defined $opts{'target'} || ! -e $opts{'xam'});
  pod2usage(-message => "Option 't|type' Type of r|target_file provided [gff3|bed].", -verbose => 0) if(!defined $opts{'type'} || ! grep ($opts{'type'},PCAP::Bam::Coverage::target_types()) );

  return \%opts;
}

__END__

=head1 NAME

xam_coverage_bins.pl - Script takes targets and a xam file and calculates coverage. Outputs it to a JSON string

=head1 SYNOPSIS

xam_coverage_bins.pl [options]

  Required parameters:
    -xam_file              -f    bam|cram file to check coverage.
    -target_file           -r    bed|gff3 file of targets.
    -output_file           -o    file to write JSON string output of coverage
    -type                  -t    Type of target file provided [bed|gff3]

  Other:
    -version               -v   Print version and exit.
    -help                  -h   Brief help message.
    -man                   -m   Full documentation.

=head1 OPTIONS

=over 8

=item B<-type>

Type of target file passed [bed|gff]

=item B<-xam_files>

bam|cram file to check coverage.

=item B<-target_file>

bed|gff3 file of targets.

=item B<-help>

Prints the help for this script

=item B<-man>

Prints the man page for this script

=back

=head1 DESCRIPTION

B<xam_coverage_bins.pl>

=cut

