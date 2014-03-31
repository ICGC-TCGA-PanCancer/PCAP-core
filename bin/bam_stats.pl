#!/usr/bin/env perl

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
use warnings FATAL => 'all';
use autodie qw(:all);

use Getopt::Long;
use File::Spec;
use Try::Tiny;
use Pod::Usage qw(pod2usage);

use PCAP::Bam::Stats;

{
  my $opts = setup();
  my $output;
  my $plots_dir = $opts->{'plots'};
  my $out_location;

  if($opts->{'output'}){
    open($output, '>', $opts->{'output'}) or die "Cannot open |".$opts->{'output'}."| for reading: $!";
    $out_location = $opts->{'output'};
  }else{
    $output = \*STDOUT;
    $out_location = 'STDOUT';
  }

  try{
    my $stats = new PCAP::Bam::Stats(-path => $opts->{'input'}, -qscoring => defined $plots_dir);
    $stats->bas($output);
    $stats->fqplots($plots_dir) if($plots_dir);
  }catch{
    die 'Reading: |'.$opts->{'input'}."| Writing to: |$out_location| Error: $_";
  }finally{
    close $output or die "Unable to close |".$opts->{'output'}."|: $!" if($opts->{'output'});
  }
}

sub setup{
  my %opts;
  my @random_args;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'i|input=s' => \$opts{'input'},
              'o|output=s' => \$opts{'output'},
              'p|plots=s' => \$opts{'plots'},
              '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  my $version = PCAP::Bam::Stats->VERSION;

  if(defined $opts{'v'}){
    print "Version: $version\n";
    exit;
  }

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  pod2usage(-message  => "\nERROR: unrecognised commandline arguments: ".join(', ',@random_args).".\n", -verbose => 1,  -output => \*STDERR) if(scalar @random_args) ;
  pod2usage(-message  => "\nERROR: i|input must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'input'});
  pod2usage(-message  => "\nERROR: i|input |".$opts{'input'}."| must be a valid file.\n", -verbose => 1,  -output => \*STDERR) unless(-f $opts{'input'});

  return \%opts;
}

__END__

=head1 NAME

bam_stats.pl - Generates a file containing read statistics for a given bam file.

=head1 SYNOPSIS

bam_stats.pl [options] [file...]

  Required parameters:
    -input    -i   File path to read in.
    -output   -o   File path to output. Defaults to STDOUT.

  Optional parameters:
    -plots    -p   Folder to contain quality score plots.

  Other:
    -help     -h   Brief help message.
    -man      -m   Full documentation.
    -version  -v   Prints the version number.

    bam_stats.pl -i my.bam -o my.bam.bas
    cat my.bam | bam_stats.pl > my.bam.bas

=head1 OPTIONS

=over 8

=item B<-input>

File path to read. Accepts only .bam files.

=item B<-output>

File path to output data. If this option is omitted the script will attempt to write to
STDOUT. The

=item B<-plots>

Directory to place quality plot images. If omitted no information about base qualities
will be collected, thus speeding up the stats collection.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-version>

Prints the version number and exits.

=back

=head1 DESCRIPTION

B<bam_stats.pl> will attempt to generate a file of statistics for a given .bam file.
The stats are are outputted as a tab-delimited rows containing a single header line.
The output consists of the following columns:

      'bam_filename' = the name of the bam file stats have been collected for.
      'sample' = the name of the sample (taken from the bam file).
      'platform' = the name of the hardware platform (taken from the bam file).
      'platform_unit' = the platform unit (i.e. lane/run) of the hardware platform (taken from the bam file).
      'library' = the library name associated with the read group.
      'readgroup' = the read group name.
      'read_length_r1' = the read length associated with read 1.
      'read_length_r2' = the read length associated with read 2.
      '#_mapped_bases' = the total number of mapped bases.
      '#_mapped_bases_r1' = the total number of mapped bases for all read 1s.
      '#_mapped_bases_r2' = the total number of mapped bases for all read 2s.
      '#_divergent_bases' = the total number of bases divergent from the reference.
      '#_divergent_bases_r1' = the total number of bases divergent from the reference for all read 1s.
      '#_divergent_bases_r2' = the total number of bases divergent from the reference for all read 2s.
      '#_total_reads' = the total number of reads.
      '#_total_reads_r1' = the total number of read 1s.
      '#_total_reads_r2' = the total number of read 2s.
      '#_mapped_reads' = the total number of unmapped reads.
      '#_mapped_reads_r1' = the total number of unmapped read 1s.
      '#_mapped_reads_r2' = the total number of unmapped read 2s.
      '#_mapped_reads_properly_paired' = the total number of properly paired reads.
      '#_gc_bases_r1' = the total number of G/C bases in read 1s.
      '#_gc_bases_r2' = the total number of G/C bases in read 2s.
      'mean_insert_size' = the mean insert size.
      'insert_size_sd' = the insert size standard deviation.
      'median_insert_size' = the median insert size.
      '#_duplicate_reads' = the total number of duplicate reads.

Any reads not linked to a ReadGroup will be combined into the group '.'
Please note that the accuracy of statistics in this group may be questionable.

=cut
