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

use Getopt::Long;
use Pod::Usage qw(pod2usage);

use Bio::DB::HTS;

const my $GFF_TYPE => 'gff';
const my $BED_TYPE => 'bed';
const my @PERMITTED_TYPES => ($GFF_TYPE,$BED_TYPE);
const my @depth_ranges => (0,1,11,21,31,41,51,101,201,501,100000000);
const my $MAX_PILEUP_DEPTH => 1_000_000;

use Data::Dumper;

my $opts = parse_options();

my $targets = parse_targets_file($opts);
my $depth_data = build_depth($opts,$targets);

if(defined($opts->{'out'})){
  my $OUT;
  open($OUT, '>', $opts->{'out'}) or croak("Error opening output file for write :$!");
    print $OUT $depth_data,"\n";
  close($OUT);
}else{
  print $depth_data,"\n";
}

sub build_depth{
  my ($opts,$targets) = @_;
  my $sam = Bio::DB::HTS->new(-bam => $opts->{'xam'});
	$sam->max_pileup_cnt($MAX_PILEUP_DEPTH);
	my $total_bases = 0;
	my %depth_bins;
  foreach my $ref_ex(@{$targets}) {
		my ($chr, $start, $end) = @{$ref_ex};
		$total_bases += $end - $start + 1;
		my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr, -start=>$start, -end=>$end);
		if(!$coverage) {
			next;
		}
		foreach my $depth(@{$coverage->coverage}) {
			if($depth_bins{$depth}) {
				$depth_bins{$depth}++;
			}
			else {
				$depth_bins{$depth} = 1;
			}
		}
	}
	my $depth_arr = build_final_bins(\%depth_bins, $total_bases);
	my $depth_data = join ',', @{$depth_arr};
	return $depth_data;
}

sub build_final_bins {
	my ($depth_bins, $total_bases) = @_;
	my %final_bins;
	my @bin_order;
	my @local_depth_ranges = @depth_ranges;
	for(my $i=0; $i<@local_depth_ranges; $i++) {
		if($i+1 < @local_depth_ranges) {
			my ($greater_equal, $less_than) = ($local_depth_ranges[$i], $local_depth_ranges[$i+1]);
			my $key = ''.$greater_equal;
			if($key ne '0') {
				if($i+2 == @local_depth_ranges) {
					$key .= '+';
				}
				elsif($greater_equal != $less_than) {
					$key .= '-'.($less_than-1);
				}
				# 0 can not be directly determined, only by doing 1 - depth[1-5]
				push @bin_order, $key;
			}

			$final_bins{$key} = 0;
			foreach my $sub_key(keys %{$depth_bins}) {
				if($key eq '0') {
					$final_bins{$key} = $depth_bins->{$sub_key};
				}
				elsif($sub_key >= $greater_equal) {
					$final_bins{$key} += $depth_bins->{$sub_key};
				}
			}
		}
	}

	my @output;
	foreach my $key(@bin_order) {
		if($key =~ m/^1\-/) {
			push @output, '0:'.(1 - sprintf('%.4f',($final_bins{$key} / $total_bases)));
		}
		push @output, $key.':'.sprintf('%.4f',($final_bins{$key} / $total_bases));
	}
	return \@output;
}

sub parse_targets_file{
  my ($opts) = @_;
  my $bait_clob;
  if($opts->{'type'} eq $GFF_TYPE){
    $bait_clob = parse_gff($opts->{'target'});
  }elsif($opts->{'type'} eq $BED_TYPE){
    $bait_clob = parse_bed($opts->{'target'});
  }

  # add a trailing \n
  $bait_clob .= "\n";
  # remove any duplicate \n
  $bait_clob =~ s/\n\n/\n/xmsg;

  my @segments;
  while(1) {
  last unless($bait_clob =~ s/([^\t]+)\t([[:digit:]]+)\t([[:digit:]]+)\n//xms);
    my ($chr, $z_start, $o_end) = ($1,$2,$3);

    croak "Start and end positions are the same, not bed format: $chr, $z_start, $o_end" if($z_start == $o_end);
    croak "Start greater than end position, not bed format: $chr, $z_start, $o_end" if($z_start > $o_end);
    push @segments, ["$chr", $z_start+1, $o_end+0]; # force data types
  }
  croak "Data remains in bait clob string after processing:\n###\n$bait_clob\n###\n" if(length $bait_clob > 0);
  return \@segments;
}

sub parse_bed{
  my ($file) =@_;
  my $FH;
  my $content;
  open($FH, '<', $file) or croak("Error trying to open file to read targets from bed: $!");
    while(<$FH>){
      next if($_ =~ m/^\s*#/); #Skip comment lines
      $content .= $_;
    }
  close($FH)  or croak("Error trying to close after reading targets from bed: $!");
  return $content;
}

sub parse_gff{
  my ($file) =@_;
  my $FH;
  my $content;
  open($FH, '<', $file) or croak("Error trying to open file to read targets from gff: $!");
    while(<$FH>){
      next if($_ =~ m/^\s*#/); #Skip comment lines
      my ($contig,$src,$type,$start,$end,undef) = split("\t",$_);
      $content .= "$contig\t".($start-1)."\t$end\n"; #Adjust start to be in bed format.
    }
  close($FH)  or croak("Error trying to close after reading targets from gff: $!");
  return $content;
}

sub parse_options{
my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'f|xam_file=s' => \$opts{'xam'},
					'r|target_file=s' => \$opts{'target'},
					'o|output_file=s' => \$opts{'out'},
					't|type=s' => \$opts{'type'},
					) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  pod2usage(-message => "Option 'f|xam_file' bam|cram file required.", -verbose => 0) if(!defined $opts{'xam'} || ! -e $opts{'xam'});
  pod2usage(-message => "Option 'r|target_file' target gff3|bed file required.", -verbose => 0) if(!defined $opts{'target'} || ! -e $opts{'xam'});
  pod2usage(-message => "Option 't|type' Type of r|target_file provided [gff3|bed].", -verbose => 0) if(!defined $opts{'type'} || ! grep ($opts{'type'},@PERMITTED_TYPES) );

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

