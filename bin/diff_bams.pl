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

use strict;
use autodie qw(:all);
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use PCAP::Cli;
use FindBin qw($Bin);
use lib "$Bin/../lib";

use Bio::DB::HTS;
use Bio::DB::HTS::AlignWrapper;

my ($file_a, $file_b, $skip_z, $count_flag_diff) = setup();

my $sam_a = Bio::DB::HTS->new(-bam => $file_a);
my $bam_a = $sam_a->hts_file;
my $bam_b = Bio::DB::HTSfile->open($file_b);

my $header_a       = $bam_a->header_read;
my $target_count_a = $header_a->n_targets;
my $target_names_a = $header_a->target_name;

my $header_b       = $bam_b->header_read;
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
my $flag_diffs = 0;
$|++;
my $last_coord = 0;
my %dup_pileup;
while(1) {
  $count++;
  $align_a = $bam_a->read1($header_a);
  $align_b = $bam_b->read1($header_b);
  if($skip_z) {
    while($align_a->qual == 0) { $align_a = $bam_a->read1($header_a); last unless(defined $align_a); }
    while($align_b->qual == 0) { $align_b = $bam_b->read1($header_b); last unless(defined $align_b); }
  }

  last if(!defined $align_a && !defined $align_b);
  die "Files have different number of records\n" if((!defined $align_a && defined $align_b) || (defined $align_a && !defined $align_b));

  die sprintf "Files differ at record $count (qname) a=%s b=%s\n", $align_a->qname, $align_b->qname if($align_a->qname ne $align_b->qname);
  if($align_a->flag ne $align_b->flag) {
    if($count_flag_diff) {
      $flag_diffs++;
      my $pos = $align_a->pos;
      if(($pos - $last_coord) >= 0) { # saves checking for different chr
        my $aw = Bio::DB::HTS::AlignWrapper->new($align_a, $sam_a);
        $dup_pileup{$aw->seq_id}{$pos} += 1;
      }
      $last_coord = $pos;
    }
    else {
      die sprintf "Files differ at record $count (flags) a=%s b=%s\n", $align_a->qname, $align_b->qname
    }
  }

  if($count % 1_000_000 == 0) {
    my $message = "Matching records: $count";
    $message .= "\t(flag mismatch: $flag_diffs)"if($count_flag_diff);
    print $message."\r";
  }
}
print "Matching records: $count\n";
if($count_flag_diff) {
  print "Flag mismatches: $flag_diffs\n";
  print "Locations of flag differences:\n";
  print "#Chr\tPos\tCount\n";
  for my $chr(sort keys %dup_pileup) {
    for my $pos(sort {$a<=>$b} keys %{$dup_pileup{$chr}}){
      print sprintf "%s\t%s\t%s\n", $chr, $pos, $dup_pileup{$chr}{$pos};
    }
  }
}

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'a|bam_a=s' => \$opts{'a'},
              'b|bam_b=s' => \$opts{'b'},
              'r|ref=s' => \$opts{'r'},
              'c|count' => \$opts{'c'},
              's|skip' => \$opts{'s'},
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  my $defined = 0;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  if($defined == 0 && @ARGV == 2) {
    $opts{'a'} = shift @ARGV;
    $opts{'b'} = shift @ARGV;
    $defined = 2;
  }
  pod2usage(-msg  => "\nERROR: Options 'a' & 'b' must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined >= 2);

  die pod2usage(-msg  => "\nERROR: bam_a and bam_b cannot be the same file\n", -verbose => 1,  -output => \*STDERR) if($opts{'a'} eq $opts{'b'});

  PCAP::Cli::file_for_reading('bam_a', $opts{'a'});
  PCAP::Cli::file_for_reading('bam_a', $opts{'b'});
  if($opts{'a'} =~ m/[.]cram/ || $opts{'b'} =~ m/[.]cram/) {
    die pod2usage(-msg  => "\nERROR: Option '-r' must be defined if any CRAM files are provided\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{'r'});
    PCAP::Cli::file_for_reading('ref', $opts{'r'});
  }

  return ($opts{'a'}, $opts{'b'}, $opts{'s'}, $opts{'c'});
}

__END__

=head1 NAME

diff_bams.pl - Compares two BAM/CRAM files irrespective of irrelevant header file differences.

=head1 SYNOPSIS

diff_bams.pl 1.bam 2.bam
  or
diff_bams.pl -a 1.bam -b 2.bam

It is rare that BAM files can be compared via MD5 due to various processing information being
incorporated into the header such as input/output paths.  This program works around this by only
comparing relevant information such as:

  1. Same number and order of @SQ headers.
  2. All reads are in same order.

You are able to compare a BAM vs it's corresponding CRAM

  Required parameters:
    -bam_a    -a    The first BAM|CRAM file.
    -bam_b    -b    The second BAM|CRAM file.

  Other:
    -ref      -r    Required for CRAM, genome.fa with co-located fai.
    -count    -c    Count flag differences
    -skipz    -s    Don't include reads with MAPQ=0 in comparison
    -help     -h    Brief help message.
    -man      -m    Full documentation.

  Run with '-m' for more detail

=head1 OPTIONS

=over 8

=item B<-bam_a>

A valid readable BAM file

=item B<-bam_b>

A valid readable BAM file (different to bam_a)

=item B<-skipz>

As reads with MAPQ zero can end up in multiple locations don't consider them
as a mismatch.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back
