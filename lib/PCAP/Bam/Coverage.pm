package PCAP::Bam::Coverage;

use strict;
use warnings FATAL=>'all';
use autodie;
use Carp;
use Const::Fast;

use Bio::DB::HTS;

const my $GFF_TYPE => 'gff3';
const my $BED_TYPE => 'bed';
const my @PERMITTED_TYPES => ($GFF_TYPE,$BED_TYPE);
const my @depth_ranges => (0,1,11,21,31,41,51,101,201,501,100000000);
const my $MAX_PILEUP_DEPTH => 1_000_000;

sub new {
  my ($class, $options) = @_;
  my $self = { };
  bless $self, $class;
  $self->_init($options);
  return $self;
}

sub _init {
  my ($self, $options) = @_;
  my $hts = Bio::DB::HTS->new(-bam => $options->{'xam'});
  $hts->max_pileup_cnt($MAX_PILEUP_DEPTH);
  $self->{'hts'} = $hts;
  $self->{'targets'} = parse_targets_file($options);
  return 1;
}

sub target_types {
  return @PERMITTED_TYPES;
}

sub build_depth{
  my $self = shift;
	my $total_bases = 0;
	my %depth_bins;
  foreach my $ref_ex(@{$self->{'targets'}}) {
		my ($chr, $start, $end) = @{$ref_ex};
		$total_bases += $end - $start + 1;
		my ($coverage) = $self->{'hts'}->features(-type=>'coverage',-seq_id=>$chr, -start=>$start, -end=>$end);
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
  my $baits;
  if($opts->{'type'} eq $GFF_TYPE){
    $baits = parse_gff($opts->{'target'});
  }elsif($opts->{'type'} eq $BED_TYPE){
    $baits = parse_bed($opts->{'target'});
  }

#use Data::Dumper;
#warn Dumper($baits);

  my @segments;
  for my $seg(@{$baits}) {
    my ($chr, $z_start, $o_end) = @{$seg};
    croak "Start and end positions are the same, not bed format: $chr, $z_start, $o_end" if($z_start == $o_end);
    croak "Start greater than end position, not bed format: $chr, $z_start, $o_end" if($z_start > $o_end);
    push @segments, ["$chr", $z_start+1, $o_end+0]; # force data types
  }

  return \@segments;
}

sub parse_bed{
  my ($file) =@_;
  my $FH;
  my @content;
  open($FH, '<', $file) or croak("Error trying to open file to read targets from bed: $!");
    while(<$FH>){
      next if($_ =~ m/^\s*#/); #Skip comment lines
      chomp $_;
      croak "File doesn't appear to be BED formatted: $_" unless($_ =~ m/^([^\t]+)\t([[:digit:]]+)\t([[:digit:]]+)/);
      push @content, [$1, $2, $3];
    }
  close($FH)  or croak("Error trying to close after reading targets from bed: $!");
  return \@content;
}

sub parse_gff{
  my ($file) =@_;
  my $FH;
  my @content;
  open($FH, '<', $file) or croak("Error trying to open file to read targets from bed: $!");
    while(<$FH>){
      next if($_ =~ m/^\s*#/); #Skip comment lines
      chomp $_;
      croak "File doesn't appear to be GFF3 formatted: $_" unless($_ =~ m/^([^\t]+)\t.+\t.+\t([[:digit:]]+)\t([[:digit:]]+)/);
      push @content, [$1, $2-1, $3];
    }
  close($FH)  or croak("Error trying to close after reading targets from bed: $!");
  return \@content;
}

1;
