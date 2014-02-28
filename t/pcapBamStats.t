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
use File::Temp qw(tempdir);
use Data::Dumper;

use Bio::DB::Sam;

const my $MODULE => 'PCAP::Bam::Stats';

my $test_data = "$Bin/../testData";
my $test_bam_file = join('/',$test_data,'Stats.bam');

## expected data... prolly needs to be const
my ($rgAnon, $rg1, $rg2) = ('.','29976','29978');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  my $obj = new_ok($MODULE => [-path => $test_bam_file]); ## this will fail without a valid sam file/object

  $obj = new_ok($MODULE => [-path => $test_bam_file, -sam => Bio::DB::Sam->new(-bam => $test_bam_file)]); ## this will fail without a valid sam file/object
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

  subtest 'read_groups' => sub {
    my $test_obj = _create_test_object();
    is_deeply($test_obj->read_groups,[$rg1, $rg2],'read_groups');
  };

  subtest 'read_group_info' => sub {
    my $test_obj = _create_test_object();
    is($test_obj->read_group_info($rg1),'@RG	ID:29976	PL:GAII	PU:5178_6	LB:PD1234a 140546_1054	PI:404	DS:short	MI:582	SM:PD1234a	PG:29976	CN:SANGER','read_group_info rg1');
    is($test_obj->read_group_info($rg2),'@RG	ID:29978	PL:GAII	PU:5085_6	LB:PD1234a 140546_1054	PI:404	DS:short	MI:582	SM:PD1234a	PG:29978	CN:SANGER','read_group_info rg2');
  };

  subtest 'calc_frac__rg' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->calc_frac_properly_paired_rg($rg1)),'0.8000','calc_frac_properly_paired_rg rg1');
    is(sprintf('%.4f',$test_obj->calc_frac_properly_paired_rg($rg2)),'0.3333','calc_frac_properly_paired_rg rg2');
    is(sprintf('%.4f',$test_obj->calc_frac_unmapped_rg($rg1)),'0.2000','calc_frac_unmapped_rg rg1');
    is(sprintf('%.4f',$test_obj->calc_frac_unmapped_rg($rg2)),'0.6667','calc_frac_unmapped_rg rg2');
    is(sprintf('%.4f',$test_obj->calc_frac_duplicate_reads_rg($rg1)),'0.4000','calc_frac_duplicate_reads_rg rg1');
    is(sprintf('%.4f',$test_obj->calc_frac_duplicate_reads_rg($rg2)),'0.0000','calc_frac_duplicate_reads_rg rg2');
  };

  subtest 'calc_frac_' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->calc_frac_properly_paired()),'0.6250','calc_frac_properly_paired');
    is(sprintf('%.4f',$test_obj->calc_frac_unmapped()),'0.3750','calc_frac_unmapped');
    is(sprintf('%.4f',$test_obj->calc_frac_duplicate_reads()),'0.2500','calc_frac_duplicate_reads');
  };

  subtest 'med_insert_size_rg' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->med_insert_size_rg($rg1),125,'med_insert_size_rg');
    is($test_obj->med_insert_size_rg($rg2),100,'med_insert_size_rg');
  };

  subtest 'med_insert_size' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->med_insert_size(),100,'med_insert_size');
  };

  subtest 'mean_insert_size_rg' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->mean_insert_size_rg($rg1)),'195.8333','mean_insert_size_rg');
    is(sprintf('%.4f',$test_obj->mean_insert_size_rg($rg2)),'100.0000','mean_insert_size_rg');
  };

  subtest 'mean_insert_size' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->mean_insert_size()),'182.1429','mean_insert_size');
  };

  subtest 'insert_sizes' => sub {
    my $test_obj = _create_test_object();
    my $expected_insert_sizes_rg1 = {
      '100' => 3,
      '150' => 1,
      '350' => 1,
      '375' => 1
    };

    my $expected_insert_sizes_rg2 = {
      '100' => 1,
    };

    is_deeply($test_obj->insert_sizes($rg1),$expected_insert_sizes_rg1,'insert_sizes rg1');
    is_deeply($test_obj->insert_sizes($rg2),$expected_insert_sizes_rg2,'insert_sizes rg2');
  };

  subtest 'properly_mapped_ratio' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->properly_mapped_ratio()),'1.6667','properly_mapped_ratio');
  };

  subtest 'read_length' => sub{
    my $test_obj = _create_test_object();

## All the reads are the same length so the output will be the same.
## Could do with adding some more gr but will change some to totals above...

    is($test_obj->read_length($rg1,1),20,'read_length rg1 read1');
    is($test_obj->read_length($rg1,2),20,'read_length rg1 read2');
    is($test_obj->read_length($rg2,1),20,'read_length rg2 read1');
    is($test_obj->read_length($rg2,2),20,'read_length rg2 read2');
  };

  subtest 'count_gc_rg' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->count_gc_rg($rg1,1),63,'count_gc_rg rg1 read1');
    is($test_obj->count_gc_rg($rg1,2),27,'count_gc_rg rg1 read2');
    is($test_obj->count_gc_rg($rg2,1),27,'count_gc_rg rg2 read1');
    is($test_obj->count_gc_rg($rg2,2),27,'count_gc_rg rg2 read2');
  };

  subtest 'mean_gc_rg' => sub {
    my $test_obj = _create_test_object();

## All the reads are the same length so the output will be the same.
## Could do with adding some more gr but will change some to totals above...
## This has already caused one bug to slip through...

    is($test_obj->mean_gc_rg($rg1,1),0.45,'mean_gc_rg rg1 read1');
    is($test_obj->mean_gc_rg($rg1,2),0.45,'mean_gc_rg rg1 read2');
    is($test_obj->mean_gc_rg($rg2,1),0.45,'mean_gc_rg rg2 read1');
    is($test_obj->mean_gc_rg($rg2,2),0.45,'mean_gc_rg rg2 read2');
  };


  subtest '_count__rg' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->count_total_reads_rg($rg1,1),7,'count_total_reads_rg rg1 read1');
    is($test_obj->count_total_reads_rg($rg1,2),3,'count_total_reads_rg rg1 read2');
    is($test_obj->count_total_reads_rg($rg2,1),3,'count_total_reads_rg rg2 read1');
    is($test_obj->count_total_reads_rg($rg2,2),3,'count_total_reads_rg rg2 read2');

    is($test_obj->count_properly_paired_rg($rg1,1),6,'count_properly_paired_rg rg1 read1');
    is($test_obj->count_properly_paired_rg($rg1,2),2,'count_properly_paired_rg rg1 read2');
    is($test_obj->count_properly_paired_rg($rg2,1),1,'count_properly_paired_rg rg2 read1');
    is($test_obj->count_properly_paired_rg($rg2,2),1,'count_properly_paired_rg rg2 read2');

    is($test_obj->count_unmapped_rg($rg1,1),1,'count_unmapped_rg rg1 read1');
    is($test_obj->count_unmapped_rg($rg1,2),1,'count_unmapped_rg rg1 read2');
    is($test_obj->count_unmapped_rg($rg2,1),2,'count_unmapped_rg rg2 read1');
    is($test_obj->count_unmapped_rg($rg2,2),2,'count_unmapped_rg rg2 read2');

    is($test_obj->count_duplicate_reads_rg($rg1,1),4,'count_duplicate_reads_rg rg1 read1');
    is($test_obj->count_duplicate_reads_rg($rg1,2),0,'count_duplicate_reads_rg rg1 read2');
    is($test_obj->count_duplicate_reads_rg($rg2,1),0,'count_duplicate_reads_rg rg2 read1');
    is($test_obj->count_duplicate_reads_rg($rg2,2),0,'count_duplicate_reads_rg rg2 read2');
  };

  subtest '_count_' => sub {
    my $test_obj = _create_test_object();

    is($test_obj->count_total_reads(1),10,'count_total_reads read1');
    is($test_obj->count_total_reads(2),6,'count_total_reads read2');

    is($test_obj->count_properly_paired(1),7,'count_properly_paired read1');
    is($test_obj->count_properly_paired(2),3,'count_properly_paired read2');

    is($test_obj->count_unmapped(1),3,'count_unmapped read1');
    is($test_obj->count_unmapped(2),3,'count_unmapped read2');

    is($test_obj->count_duplicate_reads(1),4,'count_duplicate_reads read1');
    is($test_obj->count_duplicate_reads(2),0,'count_duplicate_reads read2');
  };

  subtest 'cal_insert_size_sd_rg' => sub {
    my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->cal_insert_size_sd_rg($rg1)),'119.3879','cal_insert_size_sd_rg rg1');
    is(sprintf('%.4f',$test_obj->cal_insert_size_sd_rg($rg2)),'0.0000','cal_insert_size_sd_rg rg2');
  };

  subtest 'cal_insert_size_sd' => sub {
     my $test_obj = _create_test_object();

    is(sprintf('%.4f',$test_obj->cal_insert_size_sd()),'115.5069','cal_insert_size_sd');
  };

};

sub _create_test_object{
  return new $MODULE(-path=>$test_bam_file);
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

IL29_5178:2:54:17473:17010	579	1	9993	0	*	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976

IL29_5178:3:48:10131:17011	1091	1	9993	0	*	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:3:48:10131:17012	1091	1	9993	0	*	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:3:48:10131:17013	1091	1	9993	0	*	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:3:48:10131:17014	1091	1	9993	0	*	=	9993	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976

IL29_5178:5:103:3067:17015	323	1	9993	0	*	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978


## good read pairs
IL29_5178:2:54:17473:18172	67	1	9993	0	*	=	9993	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:2:54:17473:18172	131	1	9993	0	*	=	9993	-375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976

IL29_5178:3:48:10131:17417	67	1	9993	0	*	=	9993	350	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:3:48:10131:17417	131	1	9993	0	*	=	9993	350	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976

IL29_5178:5:103:3067:17734	67	1	9993	0	*	=	9993	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978
IL29_5178:5:103:3067:17734	131	1	9993	0	*	=	9993	-100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978

## one read unmapped
IL29_5178:7:35:17751:10872	69	1	9994	20	108M	=	9994	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976
IL29_5178:7:35:17751:10872	133	1	9994	0	*	=	9994	375	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29976

IL20_5085:6:61:7856:9407	69	1	9996	0	*	=	9996	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978
IL20_5085:6:61:7856:9407	133	1	9996	0	*	=	9996	100	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978

IL29_5178:8:46:3804:17877	69	1	9996	0	*	=	9996	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978
IL29_5178:8:46:3804:17877	133	1	9996	0	*	=	9996	150	CTCTTCCGATCTTTAGGGTT	;?;??>>>>F<BBDEBEEF	RG:Z:29978
