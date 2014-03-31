package PCAP::Bam::Stats;

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

use PCAP;
our $VERSION = PCAP->VERSION;

use strict;
use English qw( -no_match_vars );
use warnings FATAL=>'all';
use autodie qw( :all );
use Carp qw(croak carp);

use Const::Fast qw( const );
use Try::Tiny;
use File::Basename;

use Math::Gradient;
use List::Util qw(sum first);
use Bio::DB::Sam;
use GD::Image;

const my $PAIRED     => 1; # not needed
const my $PROPER     => 2;
const my $UNMAPPED   => 4;
const my $M_UNMAP    => 8;
const my $REVERSED   => 16;
const my $M_REVERSE  => 32;
const my $FIRST      => 64;
const my $SECOND     => 128;
const my $NON_PRI    => 256;
const my $V_FAIL     => 512;
const my $DUPLICATE  => 1024;
const my $SUPPLIMENT => 2048;

const my $INSERT_SMOOTH => 15;

const my @PAIRED_PROPERTIES => qw(proper);

1;

 sub new{
  my ($proto, @args) = @_;
  my $class = ref($proto) || $proto;
  my $self = { };
  bless $self, $class;
  $self->init(@args);
  return $self;
}

sub init{
  my($self,%args) = @_;

  my $path = $args{-path};
  my $sam  = $args{-sam};
  my $q_scoring = $args{-qscoring};

  unless ($sam && ref $sam eq 'Bio::DB::Sam'){
    $sam = Bio::DB::Sam->new(-bam => $path);
  }

  my $groups = _parse_header($sam);
  _process_reads($groups,$sam,$q_scoring);
  $self->{_file_path} = $path;
  $self->{_qualiy_scoring} = $q_scoring ? 1 : 0;
  $self->{_groups} = $groups;
}

sub _parse_header {
  my ($sam) = @_;
  my @lines = split "\n", $sam->header->text;
  my %groups;
  my $found = 0;
  for(@lines) {
    next if(index($_,'@RG') != 0);
    my ($id) = "$_\t" =~ m/ID:([^\t]+)/;
    $groups{$id}{'head'} = $_;
  }
  # for anonymous
  $groups{'.'}{'head'} = "\@RG\tID:anon\tLB:anon\tSM:anon";
  return \%groups;
}

sub _process_reads {
  my ($groups, $sam, $qualiy_scoring) = @_;
  my $bam = $sam->bam;
  while (my $a = $bam->read1) {
    my $flag = $a->flag;
    next if($flag & $NON_PRI); # skip secondary hits so no double counts
    next if($flag & $V_FAIL); # skip vendor fail as generally aren't considered
    next if($flag & $SUPPLIMENT); # skip supplimentary

    my $read = ($flag & $FIRST) ? 1 : 2;

    my ($rg) = ($a->get_tag_values('RG'));
    $rg ||= q{.};

    unless(exists $groups->{$rg}->{'length_'.$read}) {
      # various initialisation of elements here
      $groups->{$rg}->{'length_'.$read} = $a->l_qseq;
      $groups->{$rg}->{'fqp_'.$read} = [];
      $groups->{$rg}->{'fqp_'.$read} = [];
      $groups->{$rg}->{'inserts'} = {} if($flag & $FIRST);
    }

    $groups->{$rg}->{'count_'.$read}++;
    $groups->{$rg}->{'dup_'.$read}++ if($flag & $DUPLICATE);

    ## GC count: Total gc bases for a read group we then need to divide by read length and read count for meaning full stats
    my $qseq = $a->qseq;
    $groups->{$rg}->{'gc_'.$read} += $qseq =~ s/[GC]//gi;

    # Quality score collection for quality score plots.
    # This is really expensive therefore it is optional.
    if($qualiy_scoring){
      _add_to_qplot($groups->{$rg}->{'fqp_'.$read}, scalar $a->qscore, ($flag & $REVERSED));
    }

    if($flag & $UNMAPPED) {
      $groups->{$rg}->{'unmap_'.$read}++;
      next;
    }

    # everything after this point must require reads are mapped

    $groups->{$rg}->{'total_mapped_bases_'.$read} += ($a->calend - $a->pos);
    # Divergence calculation: Collect stats that will allow us to calculate the the number of bases that diverge from the reference.
    #                         This requires collecting the value from the NM tag and the mapped proportion of the query string.
    my ($nm) = ($a->get_tag_values('NM'));
    $groups->{$rg}->{'total_divergent_bases_'.$read} += $nm if(defined $nm);

    # Insert size can only be calculated based on reads that are on same chr
    # so it is more sensible to generate the distribution based on $PROPER-pairs.
    # only assess read 1 as size is a factor of the pair
    if(($flag & ($PROPER|$FIRST)) == ($PROPER|$FIRST)) {
      $groups->{$rg}->{'inserts'}->{abs $a->isize}++;
      $groups->{$rg}->{'proper'}++;
    }
  }

  ## Clear out empty read groups with no data
  for my $rg(keys %{$groups}) {
    my @elements = keys %{$groups->{$rg}};
    delete $groups->{$rg} if(scalar @elements == 1);
  }
}

sub _add_to_qplot {
  my ($target, $qual_ref, $reverse) = @_;
  my @quals;
  if($reverse) {
    @quals = reverse @{$qual_ref};
  }
  else {
    @quals = @{$qual_ref};
  }
  my $pos = 0;
  for(@quals) {
    $target->[$pos++]->[$_]++;
  }
}


#####################
## Calculations
#####################

sub read_groups{
  my ($self) = @_;
  return [sort {$a cmp $b} keys %{$self->{_groups}}];
}

sub read_group_info{
  my ($self, $rg) = @_;

  if($rg && exists $self->{_groups}->{$rg}){
    return $self->{_groups}->{$rg}->{'head'};
  }
}

sub read_length{
  my ($self, $rg, $read) = @_;
  return $self->{_groups}->{$rg}->{"length_$read"};
}

sub calc_frac_properly_paired_rg{
  my ($self, $rg) = @_;
  return $self->_calc_frac_property_rg('proper',$rg);
}

sub calc_frac_properly_paired{
  my ($self) = @_;
  return $self->_calc_frac_property('proper');
}

sub calc_frac_unmapped_rg{
  my ($self, $rg) = @_;
  return $self->_calc_frac_property_rg('unmap',$rg);
}

sub calc_frac_unmapped{
  my ($self) = @_;
  return $self->_calc_frac_property('unmap');
}

sub calc_frac_duplicate_reads_rg{
  my ($self, $rg) = @_;
  return $self->_calc_frac_property_rg('dup',$rg);
}

sub calc_frac_duplicate_reads{
  my ($self) = @_;
  return $self->_calc_frac_property('dup');
}

sub _calc_frac_property_rg{
  my ($self, $prop, $rg) = @_;
  my ($pp, $tt);
  my $res = 0;
  if(first {$_ eq $prop} @PAIRED_PROPERTIES) {
    $pp = ($self->{_groups}->{$rg}->{$prop} || 0);
    $tt = ($self->{_groups}->{$rg}->{'count_1'} || 0);
  }
  else {
    $pp = ($self->{_groups}->{$rg}->{$prop.'_1'} || 0) + ($self->{_groups}->{$rg}->{$prop.'_2'} || 0);
    $tt = ($self->{_groups}->{$rg}->{'count_1'} || 0) + ($self->{_groups}->{$rg}->{'count_2'} ||0);
  }

  if($tt){
    $res = $pp / $tt;
  }
  return $res;
}

sub _calc_frac_property{
  my ($self, $prop) = @_;

  my $pp = 0;
  my $tt = 0;
  my $res = 0;

  if(first {$_ eq $prop} @PAIRED_PROPERTIES) {
    foreach my $rg (keys %{$self->{_groups}}){
      $pp += ($self->{_groups}->{$rg}->{$prop} || 0);
      $tt += ($self->{_groups}->{$rg}->{'count_1'} || 0);
    }
  }
  else {
    foreach my $rg (keys %{$self->{_groups}}){
      $pp += ($self->{_groups}->{$rg}->{$prop.'_1'} || 0) + ($self->{_groups}->{$rg}->{$prop.'_2'} || 0);
      $tt += ($self->{_groups}->{$rg}->{'count_1'} || 0) + ($self->{_groups}->{$rg}->{'count_2'} || 0);
    }
  }
  if($tt){
    $res = $pp / $tt;
  }
  return $res
}

sub count_total_reads_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('count',$rg, $read);
}

sub count_total_reads{
  my ($self, $read) = @_;
  return $self->_counts('count', $read);
}

sub count_properly_paired_rg{
  my ($self, $rg) = @_;
  return $self->{_groups}->{$rg}->{'proper'} || 0;
}

sub count_properly_paired{
  my $self = shift;
  my $ret = 0;
  foreach my $rg(keys %{$self->{_groups}}){
    $ret += $self->count_properly_paired_rg($rg);
  }
  return $ret;
}

sub count_unmapped_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('unmap',$rg, $read);
}

sub count_unmapped{
  my ($self, $read) = @_;
  return $self->_counts('unmap', $read);
}

sub count_duplicate_reads_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('dup',$rg, $read);
}

sub count_duplicate_reads{
  my ($self, $read) = @_;
  return $self->_counts('dup', $read);
}

sub count_total_mapped_bases_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('total_mapped_bases',$rg, $read);
}

sub count_total_mapped_bases{
  my ($self, $read) = @_;
  return $self->_counts('total_mapped_bases', $read);
}

sub count_total_divergent_bases_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('total_divergent_bases',$rg, $read);
}

sub count_total_divergent_bases{
  my ($self, $read) = @_;
  return $self->_counts('total_divergent_bases', $read);
}

sub count_gc_rg{
  my ($self, $rg, $read) = @_;
  return $self->_counts_rg('gc',$rg, $read);
}

sub _counts_rg{
  my ($self, $prop, $rg, $read) = @_;
  croak "\$rg should be defined" unless(defined $rg);
  croak "\$rg doesn't exist" unless(exists $self->{_groups}->{$rg});
  my $ret;
  if(defined $read){
    $ret = ($self->{_groups}->{$rg}->{$prop."_$read"} || 0);
  }else{
    $ret = ($self->{_groups}->{$rg}->{$prop.'_1'} || 0) + ($self->{_groups}->{$rg}->{$prop.'_2'} || 0);
  }
  return $ret;
}

sub _counts{
  my ($self, $prop, $read) = @_;
  my $ret = 0;
  foreach my $rg(keys %{$self->{_groups}}){
    if(defined $read){
      $ret += ($self->{_groups}->{$rg}->{$prop."_$read"} || 0);
    }else{
      $ret += ($self->{_groups}->{$rg}->{$prop.'_1'} || 0);
      $ret += ($self->{_groups}->{$rg}->{$prop.'_2'} || 0);
    }
  }
  return $ret;
}

sub mean_gc_rg{
  my ($self, $rg, $read) = @_;
  croak "\$rg should be defined" unless(defined $rg);
  croak "\$rg doesn't exist" unless(exists $self->{_groups}->{$rg});

  my $ext = q{};
  $ext = "_$read" if(defined $read);
  my $prop = "gc$ext";

  croak "Property '$prop' doesn't exist for RG '$rg'" unless(exists $self->{_groups}->{$rg}->{$prop});

  if(exists $self->{_groups}->{$rg}->{$prop}){

    my $gc = $self->{_groups}->{$rg}->{$prop} || 0;
    my $read_length = $self->{_groups}->{$rg}->{'length'.$ext} || 0;
    my $read_count = $self->{_groups}->{$rg}->{'count'.$ext} || 0;

    return $gc / ($read_length * $read_count);
  }
  return undef;
}

sub insert_sizes{
  my ($self,$rg) = @_;
  if($rg && exists $self->{_groups}->{$rg}){
    return {%{$self->{_groups}->{$rg}->{'inserts'}}}; ## make a copy of the data....
  }
  return undef;
}

sub mean_insert_size_rg{
  my ($self,$rg) = @_;
  if(exists $self->{_groups}->{$rg}){
    my $pp = 0;
    my $tt = 0;

    foreach my $insert_size (keys %{$self->{_groups}->{$rg}->{'inserts'}}){
      $pp += int($insert_size) * $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
      $tt += $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
    }
    if($tt){
      return $pp / $tt;
    }else{
      return 0;
    }
  }
  return undef;
}

sub mean_insert_size{
  my ($self) = @_;

  my $pp = 0;
  my $tt = 0;
  foreach my $rg(keys %{$self->{_groups}}){
    foreach my $insert_size (keys %{$self->{_groups}->{$rg}->{'inserts'}}){
      $pp += int($insert_size) * $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
      $tt += $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
    }
  }
  if($tt){
    return $pp / $tt;
  }else{
    return 0;
  }
}

## This method returns the median insert size for the bam file. Reads with Vender Fail and Secondary Alignments are
## not included in this calculation. The method loops through all the read groups creating a single insert size bin
## hash. As we store the counts for each observed insert size, we have to implement our own median algorithm.
sub med_insert_size{
  my ($self) = @_;

  ## generate a combined insert size hash
  my $combined_insert_bin = {};
  foreach my $rg (keys %{$self->{_groups}}){
    foreach my $insert_size (keys %{$self->{_groups}->{$rg}->{'inserts'}}){
      $combined_insert_bin->{$insert_size} += $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
    }
  }
  return _med_binned_insert_sizes($combined_insert_bin);
}

## This returns the median insert size for a given read group. Reads with Vender Fail and Secondary Alignments are
## not included in this calculation. As we store the counts for each observed insert size, we have to implement our
## own median algorithm.
sub med_insert_size_rg{
  my ($self,$rg) = @_;
  if(exists $self->{_groups}->{$rg}){
    return _med_binned_insert_sizes($self->{_groups}->{$rg}->{'inserts'});
  }
  return undef;
}

sub _med_binned_insert_sizes{
  my ($insert_bin) = @_;

  ## sort the insert sizes
  my @ordered_inserts = sort {int $a <=> int $b} keys %{$insert_bin};

  my $total_conrib = 0;
  foreach my $insert_size (keys %{$insert_bin}){
    $total_conrib += $insert_bin->{$insert_size};
  }

  if( $total_conrib ){
    my $mid_point2 = int ($total_conrib / 2);
    my $mid_point = $mid_point2 + 1;

    my $running_total = 0;
    my $i = 0;

    for( $i=0;scalar(@ordered_inserts); $i++){
      $running_total += $insert_bin->{$ordered_inserts[$i]};
      last if($running_total >= $mid_point);
    }

    if($total_conrib %2 == 0 && ( $running_total - $mid_point2 >= $insert_bin->{$ordered_inserts[$i]} )){
      #warn "Thinks is even AND split between bins ";
      return ($ordered_inserts[$i] + $ordered_inserts[$i - 1]) / 2;
    }else{
      #warn "Thinks is odd or NOT split between bins";
      return $ordered_inserts[$i];
    }
  }
  return undef;
}

sub cal_insert_size_sd_rg{
  my ($self,$rg) = @_;

  if(exists $self->{_groups}->{$rg}){
    my $mean_insert_size = $self->mean_insert_size_rg($rg);
    my $pp = 0;
    my $tt = 0;

    foreach my $insert_size (keys %{$self->{_groups}->{$rg}->{'inserts'}}){
      my $count = $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
      my $dif = int($insert_size) - $mean_insert_size;
      $pp += (( $dif * $dif ) * $count);
      $tt += $count;
    }

    if($tt){
      my $variance = $pp / $tt;
      return sqrt abs $variance;
    }else{
      return 0;
    }
  }
  return undef;
}

sub cal_insert_size_sd{
  my ($self) = @_;
  my $mean_insert_size = $self->mean_insert_size();
  my $pp = 0;
  my $tt = 0;

  foreach my $rg (keys %{$self->{_groups}}){
    foreach my $insert_size (keys %{$self->{_groups}->{$rg}->{'inserts'}}){
      my $count = $self->{_groups}->{$rg}->{'inserts'}->{$insert_size};
      my $dif = int($insert_size) - $mean_insert_size;
      $pp += (( $dif * $dif ) * $count);
      $tt += $count;
    }
  }

  if($tt){
    my $variance = $pp / $tt;
    return sqrt abs $variance;
  }else{
    return 0;
  }
}

sub properly_mapped_ratio{
  my ($self) = @_;

  my $pp = 0;
  my $um = 0;
  my $tt = 0;
  my $res = 0;

  foreach my $rg (keys %{$self->{_groups}}){
    $pp += ($self->{_groups}->{$rg}->{'proper'} || 0);
    $tt += ($self->{_groups}->{$rg}->{'count_1'} || 0);
  }

#  croak("Total read count does not match combined unmapped/properly paired count. total: $tt, combined:".($pp+$um)) if $tt != ($pp+$um);

  if($tt){
    $res = $pp / $tt;
  }
  return $res;
}


#####
# PBQ COLOUR PLOT IMAGE FUNCTIONS
#####

sub fqplots {
  my ($self,$output_dir_path) = @_;
  if($self->{_qualiy_scoring}){
    my $groups = $self->{_groups};
    for my $rg(keys %{$groups}) {
      for my $read(1..2) {
        my $plot_vals = $groups->{$rg}->{'fqp_'.$read};
        down_pop_quals($plot_vals); # adds mem bloat as undef values are all filled
        fastq2image($output_dir_path, $plot_vals, $rg, $read, $groups->{$rg}->{'length_'.$read}, $groups->{$rg}->{'count_'.$read});
        delete $groups->{$rg}->{'fqp_'.$read};
      }
    }
  }
}

sub down_pop_quals {
  my ($plot_vals) = @_;
  my $max_val = (scalar @{$plot_vals}) - 1;
  for my $pos(0..$max_val) {
    my $total_qual_count = 0;
    my $raw_qual_pos = 50;
    while(1) {
      $plot_vals->[$pos]->[$raw_qual_pos] = 0 unless(defined $plot_vals->[$pos]->[$raw_qual_pos]);
      $total_qual_count += $plot_vals->[$pos]->[$raw_qual_pos];
      $plot_vals->[$pos]->[$raw_qual_pos] = $total_qual_count;
      last if($raw_qual_pos-- == 0);
    }
  }
}

sub fastq2image {
  my ($out_dir, $qual_data, $rg, $end, $num_cycles, $read_count) = @_;

  my $out_file = sprintf "$out_dir/%s_%s.png", (($rg ne q{.}) ? $rg : 'anon'), $end;

  my $read_pct = $read_count / 100; # way more efficient than ((value / read_count) * 100) in guts

  my $read ;
  my $shift = 5;

  my $width = $num_cycles * $shift + 80;

  my @quals;
  if($end == 1) {
    $read = 'forward';
     @quals = @{$qual_data};
  }
  else {
    $read = 'reverse';
    @quals = reverse @{$qual_data};
  }

  my $height = 50 * $shift + 55;

  my $im = GD::Image->new($width, $height);
  if (!$im) {
    croak q[Failed to create an image object];
  }

  my $colours = {};
  my $white = $im->colorAllocate(255,255,255);
  $im->transparent($white);

  my $black = $im->colorAllocate(0,0,0);
  my @hot_spots = ([0, 0, 0] ,[ 255, 0, 0 ], [ 255, 255, 0 ] , [0, 0, 255], [ 0, 255, 0 ]);
  my @gradient = Math::Gradient::multi_array_gradient(101, @hot_spots);
  my @colour_table = ();

  my $count = 0;
  foreach my $g (@gradient) {
    $colours->{$count} = $im->colorAllocate($g->[0], $g->[1], $g->[2]);
    push @colour_table, $colours->{$count};
    $count++;
  }

  my ($y1,$x1,$x2,$y2) = (0,40,0,0);

  my $cycle_count = 0;
  while ($cycle_count++ < $num_cycles) {
    $y1 = 10;
    $x2 = $x1 + $shift;
    my $values = shift @quals;
    my $quality = 50;
    my $value;
    while (scalar @{$values}) {
      if ($cycle_count == 1 && ($quality == 1 || ($quality > 1 && $quality % 5 == 0))) {
        $im->string(GD::gdSmallFont, 25, $y1-5, $quality, $black);
      }
      $y2 = $y1 + $shift;
      $im->filledRectangle($x1,$y1,$x2,$y2, $colours->{int ((pop @{$values}) / $read_pct)});
      $y1 = $y2;
      $quality -= 1;
    }

    if ($cycle_count == 1 || $cycle_count % 5 == 0) {
      $im->string(GD::gdSmallFont, $x1, $y1+5, $cycle_count, $black);
    }
    $x1 = $x2;
  }

  my $xaxis_label = q[Cycle number];
  if ($read) {
    $xaxis_label .= qq[ ($read read)];
  }

  my $start_xaxis_label = (int $num_cycles*$shift/2) - 40;
  if ($start_xaxis_label < 0) { $start_xaxis_label = 0; }
  $im->string(GD::gdSmallFont, $start_xaxis_label, $y1+20, $xaxis_label, $black);

  $im->stringUp(GD::gdSmallFont, 5, $height/2, q[Quality], $black);

  open my $PNG, '>', $out_file;
  binmode($PNG);
  print $PNG $im->png;
  close $PNG;
}


#####
# Output functions
#####

sub bas{
	my($self,$output_fh) = @_;

	croak 'Undefined output_path argument' unless $output_fh;

  my $file_name = basename($self->{_file_path});

  try{
    my $header = join("\t",
      'bam_filename',
#      'md5',
#      'study',
      'sample',
      'platform',
      'platform_unit',
      'library',
      'readgroup',
      'read_length_r1',
      'read_length_r2',
#      '#_total_bases',
      '#_mapped_bases',
      '#_mapped_bases_r1',
      '#_mapped_bases_r2',
      '#_divergent_bases',
      '#_divergent_bases_r1',
      '#_divergent_bases_r2',
      '#_total_reads',
      '#_total_reads_r1',
      '#_total_reads_r2',
      '#_mapped_reads',
      '#_mapped_reads_r1',
      '#_mapped_reads_r2',
      '#_mapped_reads_properly_paired',
      '#_gc_bases_r1',
      '#_gc_bases_r2',
#      '%_of_mismatched_bases',
#      'average_quality_of_mapped_bases',
      'mean_insert_size',
      'insert_size_sd',
      'median_insert_size',
#      'insert_size_median_absolute_deviation',
      '#_duplicate_reads',
#      '#_duplicate_bases'
    )."\n";

    print $output_fh $header or croak "Unable to write header line";

    foreach my $rg (@{$self->read_groups}){

      my $rg_info = $self->read_group_info($rg);

      my ($sample)   = $rg_info =~ /\tSM:([^\t]+)/;
      my ($platform) = $rg_info =~ /\tPL:([^\t]+)/;
      my ($platform_unit) = $rg_info =~ /\tPU:([^\t]+)/;
      my ($library)  = $rg_info =~ /\tLB:([^\t]+)/;

      # all header items are potentially blank
      $sample ||= q{.};
      $platform ||= q{.};
      $platform_unit ||= q{.};
      $library ||= q{.};

      my $unmapped_reads = $self->count_unmapped_rg($rg);
      my $unmapped_reads_r1 = $self->count_unmapped_rg($rg,1);
      my $unmapped_reads_r2 = $self->count_unmapped_rg($rg,2);

      my $total_reads = $self->count_total_reads_rg($rg);
      my $total_reads_r1 = $self->count_total_reads_rg($rg,1);
      my $total_reads_r2 = $self->count_total_reads_rg($rg,2);

      my $mapped_reads = $total_reads - $unmapped_reads;

      my $gc_rg1 = $self->count_gc_rg($rg,1);
      my $gc_rg2 = $self->count_gc_rg($rg,2);

      my ($mapped_reads_r1, $mapped_reads_r2, $proper_pairs) = (0,0,0);
      my ($mapped_bases, $mapped_bases_r1, $mapped_bases_r2) = (0,0,0);
      my ($divergent_bases, $divergent_bases_r1, $divergent_bases_r2) = (0,0,0);
      my ($mean_insert_size, $insert_size_sd, $median_insert_size, $dup_reads) = qw(. . . .);

      if($mapped_reads > 0) {
        # only need to count from read 1 for pairs.
        $proper_pairs = $self->count_properly_paired_rg($rg);

        $mapped_reads_r1 = $total_reads_r1 - $unmapped_reads_r1;
        $mapped_reads_r2 = $total_reads_r2 - $unmapped_reads_r2;

        $mapped_bases = $self->count_total_mapped_bases_rg($rg);
        $mapped_bases_r1 = $self->count_total_mapped_bases_rg($rg,1);
        $mapped_bases_r2 = $self->count_total_mapped_bases_rg($rg,2);

        $divergent_bases = $self->count_total_divergent_bases_rg($rg);
        if($divergent_bases > 0) {
          $divergent_bases_r1 = $self->count_total_divergent_bases_rg($rg,1);
          $divergent_bases_r2 = $self->count_total_divergent_bases_rg($rg,2);
        }
        $mean_insert_size = sprintf('%.3f',$self->mean_insert_size_rg($rg));
        $insert_size_sd = sprintf('%.3f',$self->cal_insert_size_sd_rg($rg));
        $median_insert_size = sprintf('%.3f',$self->med_insert_size_rg($rg));
        $dup_reads = $self->count_duplicate_reads_rg($rg);
      }

      my $read_length_1 = $self->read_length($rg,1);
      my $read_length_2 = $self->read_length($rg,2);

      print $output_fh join("\t",
        $file_name,
#      'md5',
#      'study',
        $sample,
        $platform,
        $platform_unit,
        $library,
        $rg,
        $read_length_1,
        $read_length_2,
#      '#_total_bases',
        $mapped_bases,
        $mapped_bases_r1,
        $mapped_bases_r2,
        $divergent_bases,
        $divergent_bases_r1,
        $divergent_bases_r2,
        $total_reads,
        $total_reads_r2,
        $total_reads_r2,
        $mapped_reads,
        $mapped_reads_r1,
        $mapped_reads_r2,
        $proper_pairs,
        $gc_rg1,
        $gc_rg2,
#      '%_of_mismatched_bases',
#      'average_quality_of_mapped_bases',
        $mean_insert_size,
        $insert_size_sd,
        $median_insert_size,
#      'insert_size_median_absolute_deviation',
        $dup_reads,
#      '#_duplicate_bases'
      )."\n";
    }

  }catch{
    croak("Problem writing output: $_");
  };
}

