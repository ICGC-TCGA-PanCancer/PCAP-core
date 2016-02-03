use strict;
use Test::More;
use Test::Fatal;
use Test::Warn;
use File::Spec;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);


const my $MODULE => 'PCAP::Bam';
const my $COMMENT_COUNT => 2;
const my $EXPECTED_RG_PL => 'HiSeq';
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
my $test_data = "$Bin/data";

# if you need to change the header.bam
# edit header.sam then:
# samcat -b header.sam > header.bam
# bamaddtrailer header.bam
## these commands are part of samcat project, queried if new samtools rc2 can handled this

my $test_bam = File::Spec->catfile($test_data, 'header.bam');
my $multi_bam = File::Spec->catfile($test_data, 'multi_sample.bam');
my $paired_bam = File::Spec->catfile($test_data, 'paired.bam');
my $unpaired_bam = File::Spec->catfile($test_data, 'unpaired.bam');
my $no_rg_bam = File::Spec->catfile($test_data, 'no_readgroups.bam');
my $md5_bam = File::Spec->catfile($test_data, 'md5.bam');
my $md5_bam_md5 = File::Spec->catfile($test_data, 'md5.bam.md5');


subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  my $obj = new_ok($MODULE => [$md5_bam]);
  is($obj->{'md5'}, $md5_bam_md5, 'md5 file for bam found');
};

subtest 'co line checks' => sub {
  my $obj = new_ok($MODULE => [$test_bam]);
  my @comments = @{$obj->comments};
  is(scalar @comments, $COMMENT_COUNT, 'Expected number of comment lines');
  unlike($comments[0], qr/\@CO/, 'Leading tag removed');
};

subtest 'rg line checks' => sub {
  my ($rg_line, $bio_db_sam) = PCAP::Bam::rg_line_for_output($test_bam);
  is($rg_line, $EXPECTED_RGLINE, 'Retreived single RG line');
  like(exception{ PCAP::Bam::rg_line_for_output($multi_bam) }
      , qr/BAM file appears to contain data for multiple readgroups, not supported unless 'existing_rgid' is found:/m
      , 'Fail when multiple readgroups in BAM');

  is((PCAP::Bam::rg_line_for_output($multi_bam, undef, undef, 1))[0], $EXPECTED_RGLINE, 'Correct RG from multiple RG header');

  my $obj = new_ok($MODULE => [$test_bam]);
  is($obj->single_rg_value('PL'), $EXPECTED_RG_PL, 'Tag successfully retrieved');
  $obj = new_ok($MODULE => [$multi_bam]);
  like(exception{$obj->single_rg_value('PL')}
      , qr/ERROR: This BAM includes multiple readgroups/m
      , 'single_rg_value not suitable for multiple RGs');
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

  is_deeply($obj->read_group_info(), $EXPECTED_MULTI_RG, 'Successful without tag presence checks');

  like(exception { $obj->read_group_info(['XX']) }
      , qr/ not found in RG of /m
      , 'Fail when required tag not found');

  $obj = new_ok($MODULE => [$no_rg_bam] );
   like(exception { $obj->read_group_info() }
      , qr/ERROR: This BAM has no readgroups:/m
      , 'Fail when no read groups found in header');

  warning_like {($sample, $bio_db_sam) = PCAP::Bam::sample_name($no_rg_bam);} qr/WARN: Failed to find samplename in RG headers of /, "Warn when no sample found";

  like(exception{ PCAP::Bam::sample_name($no_rg_bam, 1) }
      , qr/ERROR: Failed to find samplename in RG headers of /m
      , 'Die when no sample found (and die flag set)');
};

done_testing();
