#!/usr/bin/perl

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


BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use PCAP::Cli;

use Bio::DB::Sam;

my ($file_a, $file_b, $skip_z) = setup();

my $bam_a = Bio::DB::Bam->open($file_a);
my $bam_b = Bio::DB::Bam->open($file_b);

my $header_a       = $bam_a->header;
my $target_count_a = $header_a->n_targets;
my $target_names_a = $header_a->target_name;

my $header_b       = $bam_b->header;
my $target_count_b = $header_b->n_targets;
my $target_names_b = $header_b->target_name;

unless($target_count_a == $target_count_b) {
  die "Reference sequence count is different\n";
}
print "Reference sequence count passed\n";

for(0..($target_count_a-1)) {
  die "Reference sequences in different order\n" if($target_names_a->[$_] ne $target_names_b->[$_]);
}
print "Reference sequence order passed\n";

my ($align_a, $align_b, $count);
$|++;
while(1) {
  $count++;
  $align_a = $bam_a->read1;
  $align_b = $bam_b->read1;
  if($skip_z) {
    while($align_a->qual == 0) { $align_a = $bam_a->read1; last unless(defined $align_a); }
    while($align_b->qual == 0) { $align_b = $bam_b->read1; last unless(defined $align_b); }
  }


  last if(!defined $align_a && !defined $align_b);
  die "Files have different number of records\n" if((!defined $align_a && defined $align_b) || (defined $align_a && !defined $align_b));

  die sprintf "Files differ at record $count (qname) a=%s b=%s\n", $align_a->qname, $align_b->qname if($align_a->qname ne $align_b->qname);
  die sprintf "Files differ at record $count (flags) a=%s b=%s\n", $align_a->qname, $align_b->qname if($align_a->flag ne $align_b->flag);

  print "Matching records: $count\r" if($count % 100_000 == 0);
}
print "Matching records: $count\n";
print "Files match\n";

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'a|bam_a=s' => \$opts{'a'},
              'b|bam_b=s' => \$opts{'b'},
              's|skip' => \$opts{'s'},
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  my $defined = 0;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  if(@ARGV == 2) {
    $opts{'a'} = shift @ARGV;
    $opts{'b'} = shift @ARGV;
    $defined = 2;
  }
  pod2usage(-msg  => "\nERROR: Options 'a' & 'b' must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined == 2);

  die pod2usage(-msg  => "\nERROR: bam_a and bam_b cannot be the same file\n", -verbose => 1,  -output => \*STDERR) if($opts{'a'} eq $opts{'b'});

  PCAP::Cli::file_for_reading('bam_a', $opts{'a'});
  PCAP::Cli::file_for_reading('bam_a', $opts{'b'});

  return ($opts{'a'}, $opts{'b'}, $opts{'s'});
}

__END__

=head1 NAME

diff_bams.pl - Compares two BAM files irrespective of irrelevant header file differences.

=head1 SYNOPSIS

diff_bams.pl 1.bam 2.bam
  or
diff_bams.pl -a 1.bam -b 2.bam

It is rare that BAM files can be compared via MD5 due to various processing information being
incorporated into the header such as input/output paths.  This program works around this by only
comparing relevant information such as:

  1. Same number and order of @SQ headers.
  2. All reads are in same order.

  Required parameters:
    -bam_a    -a    The first BAM file.
    -bam_b    -b    The second BAM file.

  Other:
    -help     -h    Brief help message.
    -man      -m    Full documentation.

  Run with '-m' for more detail

=head1 OPTIONS

=over 8

=item B<-bam_a>

A valid readable BAM file

=item B<-bam_b>

A valid readable BAM file (different to bam_a)

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back
