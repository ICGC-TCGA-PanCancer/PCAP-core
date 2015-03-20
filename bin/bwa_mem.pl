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
use warnings FATAL => 'all';
use autodie qw(:all);
use FindBin qw($Bin);
use lib "$Bin/../lib";

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy qw(copy move);

use PCAP::Cli;
use PCAP::Bam;
use PCAP::Bwa;
use PCAP::Bwa::Meta;
use PCAP::Threaded;
use version;

const my @VALID_PROCESS => qw(setup split bwamem mark stats);
const my %INDEX_FACTOR => ( 'setup' => 1,
                            'split' => -1,
                            'bwamem' => -1,
                            'mark'   => 1,
                            'stats'  => 1,);

{
  my $options = setup();

 	my $threads = PCAP::Threaded->new($options->{'threads'});
	&PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register processes
	$threads->add_function('split', \&PCAP::Bwa::split_in);
	$threads->add_function('bwamem', \&PCAP::Bwa::bwa_mem, exists $options->{'index'} ? 1 : &PCAP::Bwa::bwa_mem_max_cores);

  PCAP::Bwa::mem_setup($options) if(!exists $options->{'process'} || $options->{'process'} eq 'setup');

	$threads->run($options->{'max_split'}, 'split', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'split');

  if(!exists $options->{'process'} || $options->{'process'} eq 'bwamem') {
    $options->{'max_index'} = PCAP::Bwa::mem_mapmax($options);
    $threads->run($options->{'max_index'}, 'bwamem', $options);
  }

  PCAP::Bam::merge_and_mark_dup($options, File::Spec->catdir($options->{'tmp'}, 'sorted')) if(!exists $options->{'process'} || $options->{'process'} eq 'mark');
  if(!exists $options->{'process'} || $options->{'process'} eq 'stats') {
    PCAP::Bam::bam_stats($options);
    &cleanup($options);
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'j|jobs' => \$opts{'jobs'},
              't|threads=i' => \$opts{'threads'},
              'r|reference=s' => \$opts{'reference'},
              'o|outdir=s' => \$opts{'outdir'},
              's|sample=s' => \$opts{'sample'},
              'c|charts' => \$opts{'charts'},
              'f|fragment=s' => \$opts{'fragment'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  if(defined $opts{'v'}) {
    print PCAP->VERSION,"\n";
    exit 0;
  }

  my $version = PCAP::Bwa::bwa_version();
  die "bwa mem can only be used with bwa version 0.7+, the version found in path is: $version\n" unless(version->parse($version) >= version->parse('0.7.0'));

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined);

  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});

  # now safe to apply defaults
  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpMap_'.$opts{'sample'});
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;
  $opts{'raw_files'} = \@ARGV;

  my $max_split = PCAP::Bwa::mem_prepare(\%opts);

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    my $max_index = $INDEX_FACTOR{$opts{'process'}};
    if($max_index == -1) {
      $max_index = $max_split                     if($opts{'process'} eq 'split');
      $max_index = PCAP::Bwa::mem_mapmax(\%opts)  if($opts{'process'} eq 'bwamem');
    }
    $opts{'max_index'} = $max_index;
    if(exists $opts{'index'}) {
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, \@ARGV, $max_index);
    }
    elsif(defined $opts{'jobs'}) {
      # tell me the max processes required for this step
      print "Requires: $max_index\n";
      exit 0;
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  return \%opts;
}

__END__

=head1 NAME

bwa_mem.pl - Align a set of lanes to specified reference with single command.

=head1 SYNOPSIS

bwa_aln.pl [options] [file(s)...]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -reference -r   Path to reference genome file *.fa[.gz]
    -sample    -s   Sample name to be applied to output file.
    -threads   -t   Number of threads to use. [1]

  Optional parameters:
    -fragment  -f   Split input into fragements of X million repairs
    -charts    -c   Generate perl-base-quality plots (+30% runtime on stats step)

  Targeted processing:
    -process   -p   Only process this step then exit, optionally set -index
                      bwamem - only applicable if input is bam
                        mark - Run duplicate marking (-index N/A)
                       stats - Generates the *.bas file for the final BAM.

    -index     -i   Optionally restrict '-p' to single job
                      bwamem - 1..<lane_count>

  Other:
    -jobs      -j   For a parallel step report the number of jobs required
    -help      -h   Brief help message.
    -man       -m   Full documentation.

  File list can be full file names or wildcard, e.g.
    bwa_mem.pl -t 16 -r some/genome.fa.gz -o myout -s sample input/*.bam

  Run with '-m' for possible input file types.

=head1 OPTIONS

=over 8

=item B<-outdir>

Directory to write output to.  During processing a temp folder will be generated in this area,
should the process fail B<only delete this if> you are unable to resume the process.

Final output files include: <SAMPLE>.bam, <SAMPLE>.bam.bai, <SAMPLE>.md5, <SAMPLE>.met

=item B<-reference>

Path to genome.fa[.gz] file and associated indexes for BWA.

=item B<-sample>

Name to be applied to output files.  Special characters will not be magically fixed.

=item B<-threads>

Number of threads to be used in processing.

If perl is not compiled with threading some steps will not run in parallel, however much of the
script calls other tools that will still utilise this appropriately.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head2 TARGETED PROCESSING

=over 8

=item B<-process>

If you want to run the code in a more efficient manner then this allows each procesing type to be
executed in isolation.  You can restrict to a single process within the block by specifying
B<-index> as well.

=back

=head2 INPUT FILE TYPES

There are several types of file that the script is able to process.

=over 8

=item f[ast]q

A standard uncompressed fastq file.  Requires a pair of inputs with standard suffix of '_1' and '_2'
immediately prior to '.f[ast]q' or an interleaved f[ast]q file where read 1 and 2 are adjacent
in the file.


=item f[ast]q.gz

As *.f[ast]q but compressed with gzip.

=item bam

A list of single lane BAM files, RG line is transfered to aligned files.

=back

=head1 DESCRIPTION

B<bwa_mem.pl> will attempt to run all mapping steps for BWA-mem, as well as subsequent merging
and duplicate marking automatically.

=cut

