#!/usr/bin/perl

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

use Cwd qw(abs_path);
use strict;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use File::Basename;
use Carp;
use Getopt::Long;
use Pod::Usage;

use PCAP;

use Bio::DB::BigWig 'binMean','binStdev';

my %chr_stats;
my @chr_order;

{
  my $options = option_builder();

  my $wig = Bio::DB::BigWig->new(-bigwig=>$options->{'b'});
  my @chroms = $wig->features(-type=>'summary');

  for my $c (@chroms) {
    my $seqid   = $c->seq_id;
    next if(defined $options->{'r'} && $options->{'r'} ne $seqid && 'chr'.$options->{'r'} ne $seqid);
    my $start = $c->start;

    my $stats     = $c->statistical_summary(1);
    my $bin_width = $c->length/@$stats;

    my $s = shift @{$stats};

    my $mean  = binMean($s);
    my $stdev = binStdev($s);
    my $end   = $start + $bin_width-1;

    push @chr_order, $seqid;
    $chr_stats{$seqid} = {'mean' =>  binMean($s),
                          'stdev' => binStdev($s)};
    warn sprintf "%s: mean %.2f, stdev %.2f\n", $seqid, $chr_stats{$seqid}{'mean'}, $chr_stats{$seqid}{'stdev'};
  }

  open my $OFH, '>', $options->{'o'} or die "Failed to create $options->{o}: $!\n";
  for my $chr(@chr_order) {
    my $max_val = $chr_stats{$chr}{'mean'} + ($chr_stats{$chr}{'stdev'} * $options->{'s'});
    warn sprintf "%s: Max depth permitted = %d\n", $chr, $max_val;
    my $iterator = $wig->get_seq_stream(-seq_id=> $chr);
    while (my $p = $iterator->next_seq) {
      next if($p->score <= $max_val);
      printf $OFH "%s\t%d\t%d\t%d\n", $chr, $p->start-1, $p->end, $p->score;
    }
  }
  close $OFH;
}

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'    => \$opts{'h'},
		'b|bigwig=s' => \$opts{'b'},
		'o|output=s' => \$opts{'o'},
		'r|ref=s' => \$opts{'r'},
		'd|decode=s@' => \$opts{'d'},
    's|sd=n' => \$opts{'s'},
    'v|version' => \$opts{'v'},
	);

  if(defined $opts{'v'}) {
    print PCAP->VERSION,"\n";
    exit 0;
  }

	pod2usage(0) if($opts{'h'});

	pod2usage(1) if(!$opts{'b'} || !$opts{'o'});

  croak $opts{'b'}.' was not found or is empty' if(!-e $opts{'b'} || !-s $opts{'b'});

  if($opts{'d'}) {
	  if(!$opts{'r'}) {
	    croak '-d should not be defined without -r';
	  }
	  my %decode;
	  foreach my $d_str(@{$opts{'d'}}) {
	    if($d_str =~ m/^(\d+)\:(.*)$/) {
	      my $num = $1;
	      my $chr = $2;
	      $decode{$num} = $chr;
	    }
	    else {
	      croak "Decode string of $d_str is invalid see --help";
	    }
	  }
	  if(defined $decode{$opts{'r'}}) {
	    $opts{'r'} = $decode{$opts{'r'}};
	  }
	}

	my $fn = fileparse($opts{'b'});
	$fn =~ s/\.bw$//;
	$opts{'o'} .= '/' if($opts{'o'} !~ m/\/$/);
	$opts{'o'} .= $fn;
	if($opts{'r'}) {
	  $opts{'o'} .= '.'.$opts{'r'};
	}
	$opts{'o'} .= '.bed';

	if(!$opts{'s'}) {
	  $opts{'s'} = 12;
	}

	return \%opts;
}

__END__

=head1 NAME

detectExtremeDepth.pl - Generate profile of BigWig file and identify regions outside the normal range

=head1 SYNOPSIS

  General Options (list OR project must be defined):

    --bigwig    (-b)  FILE  BigWig file path
    --output    (-o)  DIR   Folder to send output to
                             - named as input file with '.tab' extension
                             - if '-r' defined '.{val}' will prefix '.bed'

  Optional:
    --ref       (-r)  STR   Restrict to this reference (mainly for testing)
                             - without 'chr' prefix, will test with and without the 'chr' for you.
    --decode    (-d)  STR   Decode -r to chromosome names (do not include 'chr')
                             e.g. -d 23:X -d 24:Y -d 25:MT
    --sd        (-s)  INT   Number of standard deviations above mean for group to be included [12]
    --help      (-h)        This message
    --version   (-v)        Version

  Examples:
    perl ~/detectExtremeDepth.pl -o someplace -b sample.bw

=cut
