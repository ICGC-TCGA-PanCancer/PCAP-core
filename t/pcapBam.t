use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);


const my $MODULE => 'PCAP::Bam';
const my $EXPECTED_RGLINE => q{@RG\tID:1\tPL:HiSeq\tPU:1_1\tLB:SAMPLE_LIBRARY\tPI:500\tDS:short\tSM:SAMPLE_NAME\tCN:SANGER};
const my $EXPECTED_SAMPLE => 'SAMPLE_NAME';
const my $EXPECTED_SINGLE_RG => [{'CN' => 'SANGER',
                                  'ID' => '1',
                                  'PU' => '1_1',
                                  'DS' => 'short',
                                  'PL' => 'HiSeq',
                                  'LB' => 'SAMPLE_LIBRARY',
                                  'SM' => 'SAMPLE_NAME',
                                  'PI' => '500'
                                } ];
const my $EXPECTED_MULTI_RG => [ {'CN' => 'SANGER',
                                  'ID' => '1',
                                  'PU' => '1_1',
                                  'DS' => 'short',
                                  'PL' => 'HiSeq',
                                  'LB' => 'SAMPLE_LIBRARY',
                                  'SM' => 'SAMPLE_NAME',
                                  'PI' => '500'
                                },
                                 {'CN' => 'SANGER',
                                  'ID' => '2',
                                  'PU' => '1_2',
                                  'DS' => 'short',
                                  'PL' => 'HiSeq',
                                  'LB' => 'ANOTHER_LIBRARY',
                                  'SM' => 'ANOTHER_SAMPLE',
                                  'PI' => '500'
                                } ];

use FindBin qw($Bin);
my $test_data = "$Bin/../testData";

my $test_bam = File::Spec->catfile($test_data, 'header.bam');
my $multi_bam = File::Spec->catfile($test_data, 'multi_sample.bam');
my $paired_bam = File::Spec->catfile($test_data, 'paired.bam');
my $unpaired_bam = File::Spec->catfile($test_data, 'unpaired.bam');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'rg_line checks' => sub {
  my ($rg_line, $bio_db_sam) = PCAP::Bam::rg_line_for_output($test_bam);
  is($rg_line, $EXPECTED_RGLINE, 'Retreived single RG line');
  like(exception{ PCAP::Bam::rg_line_for_output($multi_bam) }
      , qr/BAM file appears to contain data for multiple readgroups, not supported:/m
      , 'Fail when multiple readgroups in BAM');
};

subtest 'paired seq checks' => sub {
  my $obj = new_ok($MODULE => [$paired_bam]);
  ok($obj->check_paired, 'Paired bam ok');
  $obj = new_ok($MODULE => [$unpaired_bam]);
  like(exception{ $obj->check_paired }
      , qr/ERROR: Input BAMs should be for paired end sequencing:/m
      , 'Fail when unpaired sequencing in BAM');
};

subtest 'header_info checks' => sub {
  my ($sample, $bio_db_sam) = PCAP::Bam::sample_name($test_bam);
  is($sample, $EXPECTED_SAMPLE, 'Retrieved sample name');
  like(exception{ PCAP::Bam::sample_name($multi_bam) }
      , qr/BAM file appears to contain data for multiple samples, not supported:/m
      , 'Fail when multiple samples in BAM');

  my $obj = new_ok($MODULE => [$test_bam]);
  is_deeply($obj->read_group_info(), $EXPECTED_SINGLE_RG, 'Expected single RG structure');

  $obj = new_ok($MODULE => [$multi_bam]);
  is_deeply($obj->read_group_info(), $EXPECTED_MULTI_RG, 'Expected multi RG structure');

  is_deeply($obj->read_group_info(['SM','ID']), $EXPECTED_MULTI_RG, 'Successful with tag presence checks');

  like(exception { $obj->read_group_info(['XX']) }
      , qr/ not found in RG of /m
      , 'Fail when required tag not found');
};

done_testing();
