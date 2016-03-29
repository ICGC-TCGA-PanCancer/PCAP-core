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
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw/ :seekable /;
use Data::Dumper;

use Bio::DB::HTS;

const my $MODULE => 'PCAP::Bam::Stats';

my $test_data = "$Bin/data";
my $test_bam_file = join('/',$test_data,'Stats.bam');
my $test_bas_file = join('/',$test_data,'Stats.bam.bas');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  my $obj = new_ok($MODULE => [-path => $test_bam_file]); ## this will fail without a valid sam file/object

  $obj = new_ok($MODULE => [-path => $test_bam_file, -sam => Bio::DB::HTS->new(-bam => $test_bam_file)]); ## this will fail without a valid sam file/object
};

subtest 'Per RG checks' => sub {
  my $tmp;
  my $test_obj = _create_test_object();

  eval{
    $tmp = File::Temp->new(TEMPLATE => '/tmp/Bam_StatsXXXXXX', SUFFIX => '.bas.tmp');
    my $fname = $tmp->filename;
    $tmp->unlink_on_destroy( 1 );
    $test_obj->bas($tmp);
    close $tmp;

    my $got = tsv_to_data($fname);
    my $expected = tsv_to_data($test_bas_file);

    is_deeply($got,$expected, "Per RG checks compared to: $test_bas_file") or BAIL_OUT "Merged results will be invalid if per-RG data is";
    1;
  };if($@){
    BAIL_OUT $@;
  }
};

subtest 'Non-object funcions' => sub {

  subtest '_add_to_qplot' => sub {
    my $target = []; ##
    my $target2 = [];
    my $qual_ref1 = [1,2,3,4,5]; ## an array of qualities
    my $qual_ref2 = [1,3,4,5,4]; ## an array of qualities
    my $reverse_y = 1;
    my $reverse_n = 0;

    my $target_after_n = [
      [undef,1],[undef,undef,1],[undef,undef,undef,1],[undef,undef,undef,undef,1],[undef,undef,undef,undef,undef,1]
    ];

    my $target_after_n2 = [
      [undef,2],[undef,undef,1,1],[undef,undef,undef,1,1],[undef,undef,undef,undef,1,1],[undef,undef,undef,undef,1,1]
    ];

    my $target_after_y = [
      [undef,undef,undef,undef,undef,1],[undef,undef,undef,undef,1],[undef,undef,undef,1],[undef,undef,1],[undef,1]
    ];
    my $target_after_y2 = [
      [undef,undef,undef,undef,1,1],[undef,undef,undef,undef,1,1],[undef,undef,undef,1,1],[undef,undef,1,1],[undef,2]
    ];

    PCAP::Bam::Stats::_add_to_qplot($target, $qual_ref1, $reverse_n);

    is_deeply($target, $target_after_n, '_add_to_qplot');
    PCAP::Bam::Stats::_add_to_qplot($target, $qual_ref2, $reverse_n);
    is_deeply($target, $target_after_n2, '_add_to_qplot');

    PCAP::Bam::Stats::_add_to_qplot($target2, $qual_ref1, $reverse_y);
    is_deeply($target2, $target_after_y, '_add_to_qplot');
    PCAP::Bam::Stats::_add_to_qplot($target2, $qual_ref2, $reverse_y);
    is_deeply($target2, $target_after_y2, '_add_to_qplot');
  };

  subtest 'down_pop_quals' => sub {
    my $plot_vals = [
      [1,2,3,4],
      [5,undef,6,7,8]
    ]; ## some vals to process

    my $plot_vals_after = [
      [10,9,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
      [26,21,21,15,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ]; ## what the array_ref will look like after processing

    PCAP::Bam::Stats::down_pop_quals($plot_vals);

    is_deeply($plot_vals, $plot_vals_after, 'down_pop_quals');
  };
};

subtest 'Object funcions' => sub {

  subtest 'calc_frac_' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->calc_frac_properly_paired()),'0.6154','calc_frac_properly_paired');
    is(sprintf('%.4f',$test_obj->calc_frac_unmapped()),'0.2727','calc_frac_unmapped');
    is(sprintf('%.4f',$test_obj->calc_frac_duplicate_reads()),'0.1818','calc_frac_duplicate_reads');
  };

  subtest 'med_insert_size' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->med_insert_size(),100,'med_insert_size');
  };

  subtest 'mean_insert_size' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->mean_insert_size()),'171.8750','mean_insert_size');
  };

  subtest 'properly_mapped_ratio' => sub {
    my $test_obj = _create_test_object();
    is(sprintf('%.4f',$test_obj->properly_mapped_ratio()),'0.6154','properly_mapped_ratio');
  };

  subtest '_count_' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->count_total_reads(1),13,'count_total_reads read1');
    is($test_obj->count_total_reads(2),9,'count_total_reads read2');

    is($test_obj->count_properly_paired(),8,'count_properly_paired');

    is($test_obj->count_unmapped(1),3,'count_unmapped read1');
    is($test_obj->count_unmapped(2),3,'count_unmapped read2');

    is($test_obj->count_duplicate_reads(1),4,'count_duplicate_reads read1');
    is($test_obj->count_duplicate_reads(2),0,'count_duplicate_reads read2');

    is($test_obj->count_total_mapped_bases(1),195,'count_total_mapped_bases read1');
    is($test_obj->count_total_mapped_bases(2),117,'count_total_mapped_bases read2');

    is($test_obj->count_total_divergent_bases(1),41,'count_total_divergent_bases read1');
    is($test_obj->count_total_divergent_bases(2),30,'count_total_divergent_bases read2');
  };

  subtest 'cal_insert_size_sd' => sub {
     my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->cal_insert_size_sd()),'111.4096','cal_insert_size_sd');
  };

};


sub _create_test_object{
  return new $MODULE(-path=>$test_bam_file);
}

sub tsv_to_data {
  my $file = shift;
  my @data;
  open (my $fh, '<', $file) or fail "Could not open $file for reading: $!";
  while(my $line = <$fh>) {
    chomp $line;
    push @data, [split /\t/, $line];
  }
  close($fh);
  return \@data;
}

done_testing();

## The data in the block below should be kept in sync with the ../testData/Stats.bam file.
## This will make it easier to se what is going in the tests.

#####reads to skip...
##579
# read paired
# read mapped in proper pair
# first in pair
# read fails platform/vendor quality checks

##1091
# read paired
# read mapped in proper pair
# first in pair
# read is PCR or optical duplicate

##323
# read paired
# read mapped in proper pair
# first in pair
# supplementary alignment

#####


##67
# read paired
# read mapped in proper pair
# first in pair

##131
# read paired
# read mapped in proper pair
# second in pair

##69
# read paired
# read unmapped
# first in pair

##133
# read paired
# read unmapped
# second in pair



__DATA__
@HD	VN:1.0	GO:none	SO:coordinate
@SQ	SN:1	LN:249250621	AS:37	SP:HUMAN
@RG	ID:29976	PL:GAII	PU:5178_6	LB:PD1234a 140546_1054	PI:404	DS:short	MI:582	SM:PD1234a	PG:29976	CN:SANGER
@RG	ID:29978	PL:GAII	PU:5085_6	LB:PD1234a 140546_1054	PI:404	DS:short	MI:582	SM:PD1234a	PG:29978	CN:SANGER
@RG	ID:29979	PL:GAII	PU:5086_6	LB:PD1234a 140546_1054	PI:404	DS:short	MI:582	SM:PD1234a	PG:29979	CN:SANGER

IL29_5178:2:54:17473:17010	579	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:1

IL29_5178:3:48:10131:17011	1091	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:2
IL29_5178:3:48:10131:17012	1091	1	9993	0	3S15M2S	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:3
IL29_5178:3:48:10131:17013	1091	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:4
IL29_5178:3:48:10131:17014	1091	1	9993	0	20M	=	9993	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:5

IL29_5178:5:103:3067:17015	323	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978	NM:i:6


## good read pairs
IL29_5178:2:54:17473:18172	67	1	9993	0	9M1I10M	=	9993	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:7
IL29_5178:2:54:17473:18172	131	1	9993	0	20M	=	9993	-375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:8

IL29_5178:3:48:10131:17417	67	1	9993	0	20M	=	9993	350	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:9
IL29_5178:3:48:10131:17417	131	1	9993	0	2S17M1S	=	9993	350	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:10

IL29_5178:5:103:3067:17734	67	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978	NM:i:11
IL29_5178:5:103:3067:17734	131	1	9993	0	20M	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978	NM:i:12

## one read unmapped
IL29_5178:7:35:17751:10872	69	1	9994	20	*	=	9994	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976
IL29_5178:7:35:17751:10872	133	1	9994	0	*	=	9994	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29976

IL20_5085:6:61:7856:9407	69	1	9996	0	*	=	9996	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978
IL20_5085:6:61:7856:9407	133	1	9996	0	*	=	9996	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978

IL29_5178:8:46:3804:17877	69	1	9996	0	*	=	9996	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978
IL29_5178:8:46:3804:17877	133	1	9996	0	*	=	9996	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29978

## Just for mapped_pairs/inter_chr_pairs
# proper pair
IL29_5086:6:35:17751:10872	67	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
IL29_5086:6:35:17751:10872	131	1	9993	0	20M	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
# non-proper
IL29_5086:6:35:17751:10873	65	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
IL29_5086:6:35:17751:10873	129	1	9993	0	20M	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
# non-proper, inter-chr
IL29_5086:6:35:17751:10874	65	1	9993	0	20M	X	9993	0	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
IL29_5086:6:35:17751:10874	129	X	9993	0	20M	1	9993	0	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEFF	RG:Z:29979
