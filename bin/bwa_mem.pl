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
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);

use PCAP::Cli;
use PCAP::Bam;
use PCAP::Bwa;
use PCAP::Bwa::Meta;
use version;

my @mod_list = keys %INC;
exit 0 if(first {$_ =~ m|^Devel/Cover| } @mod_list);

const my @VALID_PROCESS => qw(bwamem mark);
const my %INDEX_FACTOR => ( 'bwamem' => 1,
                            'mark'   => 1,);

{
  my $options = setup();
  $options->{'meta_set'} = PCAP::Bwa::Meta::files_to_meta($options->{'tmp'}, $options->{'raw_files'}, $options->{'sample'});

  if($options->{'meta_set'}->[0]->fastq && !$options->{'meta_set'}->[0]->paired_fq) {
    die "BWA aln doesn't support interleaved FASTQ\n";
  }

  my $bam_count = scalar @{$options->{'meta_set'}};
  PCAP::Bwa::bwa_mem($options) if(!exists $options->{'process'} || $options->{'process'} eq 'bwamem');
  if(!exists $options->{'process'} || $options->{'process'} eq 'mark') {
    PCAP::Bam::merge_and_mark_dup($options);
    &cleanup($options);
  }
}

sub cleanup {
  my $tmpdir = shift->{'tmp'};
  remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub setup {
  my %opts;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              't|threads=i' => \$opts{'threads'},
              'r|reference=s' => \$opts{'reference'},
              'o|outdir=s' => \$opts{'outdir'},
              's|sample=s' => \$opts{'sample'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'m'});

  my $version = PCAP::Bwa::bwa_version();
  die "bwa mem can only be used with bwa version 0.7+, the version found in path is: $version\n" unless(version->parse($version) >= version->parse('0.7.0'));

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless($defined);

  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});



  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      # theres a small assumption that @ARGV is correct but that is still checked before exec
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, \@ARGV, $INDEX_FACTOR{$opts{'process'}});
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

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
  return \%opts;
}

__END__

=head1 NAME

bwa_aln.pl - Align a set of lanes to specified reference with single command.

=head1 SYNOPSIS

bwa_aln.pl [options] [file(s)...]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -reference -r   Path to reference genome file *.fa[.gz]
    -sample    -s   Sample name to be applied to output file.
    -threads   -t   Number of threads to use. [1]

  Targeted processing:
    -process   -p   Only process this step then exit, optionally set -index
                      bwamem - only applicable if input is bam
                        mark - Run duplicate marking (-index N/A)

    -index     -i   Optionally restrict '-p' to single job
                      bwamem - 1..<lane_count>

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.

  File list can be full file names or wildcard, e.g.
    bwa_aln.pl -t 16 -r some/genome.fa.gz -o myout -s sample input/*.fq.gz

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
immediately prior to '.f[ast]q'.


=item f[ast]q.gz

As *.f[ast]q but compressed with gzip.

=item bam

A list of single lane BAM files, no information is taken from the headers.

B<This method has additional processing converted to *.fq.gz to give common start point.>

=back

=head1 DESCRIPTION

B<bwa_aln.pl> will attempt to run all mapping steps for BWA ALN and subsequent duplicate marking
automatically.

=cut

