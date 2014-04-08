use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Temp qw(tempdir);

const my $MODULE => 'PCAP::SRA';
use FindBin qw($Bin);
my $real_data = "$Bin/../share";


my $obj;
subtest 'Initialisation checks' => sub {
  local $SIG{__WARN__}=sub{};
  use_ok($MODULE);
};

subtest 'CV lookup checks' => sub {
  # horrid but best way to ensure full dataset is available for tests
  my $dir = tempdir( CLEANUP => 1 );
  `cp -r $real_data/cv_tables $dir/.`;

  is(ref PCAP::SRA::create_cv_lookups($dir), 'HASH', 'Successful load of CV terms');

  # remove one of the files
  my $tmp_pa = "$dir/cv_tables/ICGC/dcc_project_code.txt";
  unlink $tmp_pa;

  like(exception{PCAP::SRA::create_cv_lookups($dir)}
      , qr/ERROR: Unable to find controlled vocabulary file for /
      , 'Die when CV file not found');

  # re-create the missing file but empty
  `touch $tmp_pa`;

  like(exception{PCAP::SRA::create_cv_lookups($dir)}
      , qr/ERROR: Controlled vocabulary file for dcc_project_code is empty: /
      , 'Die when CV file is empty');
};

subtest 'Uncategorised checks' => sub {
  my $uuid = PCAP::SRA::uuid();
  is($uuid, lc $uuid, 'Confirm uuid has been created in lower-case');

  ok(PCAP::SRA::validate_seq_type('WGS'), 'Expected sequencing type');
  like(exception{PCAP::SRA::validate_seq_type('XXX')}
    , qr/ERROR: '.*' is not a recognised sequencing\/library type/
    , 'Die when unexpected sequencing type');
};

subtest 'Donor/sample correlation checks' => sub {
  my $sra = new_ok($MODULE);
  ok($sra->validate_info(bam_ob('tumour')), 'Check tumour validates');
  ok($sra->validate_info(bam_ob('normal')), 'Check normal validates');

  my $tum_bam = bam_ob('tumour');

  # delete a required field
  delete $tum_bam->{'info'}->{'dcc_project_code'};
  like(exception{ $sra->validate_info($tum_bam) }
    , qr/Required comment field/
    , 'Die on missing required field');

  $tum_bam = bam_ob('tumour');
  # change donor
  $tum_bam->{'info'}->{'submitter_donor_id'} = 'wibble';
  like(exception{ $sra->validate_info($tum_bam) }
    , qr/Each execution of the script should be limited to samples from the same submitter_donor_id/
    , 'Die on different Donor');

  $tum_bam = bam_ob('tumour');
  # change SM UUID
  $tum_bam->{'SM'} = 'xxxx-xxxx-xxxx-xxxx-xxxxxxxxx';
  like(exception{ $sra->validate_info($tum_bam) }
    , qr/submitter_sample_id of '.*' has multiple SM UUIDs/
    , 'Die on multiple SM UUIDs for same submitter_sample_id');

  $tum_bam = bam_ob('tumour');
  # change submitter_sample_id
  $tum_bam->{'info'}->{'submitter_sample_id'} = 'PD4115c';
  like(exception{ $sra->validate_info($tum_bam) }
    , qr/SM UUID of '.*' has multiple submitter_sample_id's/
    , q{Die on multiple submitter_sample_ids' for SM UUID});

  my $norm_bam = bam_ob('normal2');
  # try to add a different normal
  like(exception{ $sra->validate_info($norm_bam) }
    , qr/Only one normal sample can be defined for a donor/
    , q{Die when is second normal added});
};

done_testing();

sub bam_ob {
  my $type = shift;
  # these can't be constants as I need to be able to mutate them between tests
  my $bam_ob_tum =  { 'exp' => 'WTSI:5206:939',
                      'info' => {
                                  'dcc_project_code' => 'BRCA-UK',
                                  'submitter_donor_id' => 'CGP_donor_1199137',
                                  'submitter_sample_id' => 'PD4115a',
                                  'dcc_specimen_type' => 'Primary tumour - solid tissue',
                                  'use_cntl' => 'f393bb05-47c6-f3a5-e040-11ac0d48452c',
                                  'total_lanes' => 22,
                                  'submitter_specimen_id' => 'CGP_specimen_1227969',
                                  'ega_sample_accession' => 'EGAN00001000592'
                                },
                      'CN' => 'WTSI',
                      'ID' => 'WTSI6485',
                      'PU' => 'WTSI:5206_3',
                      'file' => 'PD4115a/5206_3.bam',
                      'PM' => 'Illumina Genome Analyzer II',
                      'PL' => 'ILLUMINA',
                      'run' => 'WTSI:5206',
                      'LB' => 'WGS:WTSI:939',
                      'SM' => 'f393bb08-4121-cad8-e040-11ac0d484535',
                      'type' => 'WGS',
                      'PI' => '450',
                      'md5' => 'PD4115a/5206_3.bam.md5',
                      'DT' => '2010-09-10T00:00:00.000+00:00'
                    };

  my $bam_ob_norm = { 'exp' => 'WTSI:5143:938',
                      'info' => {
                                  'dcc_project_code' => 'BRCA-UK',
                                  'submitter_donor_id' => 'CGP_donor_1199137',
                                  'submitter_sample_id' => 'PD4115b',
                                  'dcc_specimen_type' => 'Normal - blood derived',
                                  'use_cntl' => 'N/A',
                                  'total_lanes' => 21,
                                  'submitter_specimen_id' => 'CGP_specimen_1227970',
                                  'ega_sample_accession' => 'EGAN00001000593'
                                },
                      'CN' => 'WTSI',
                      'ID' => 'WTSI6066',
                      'PU' => 'WTSI:5143_6',
                      'file' => 'PD4115b/5143_6.bam',
                      'PM' => 'Illumina Genome Analyzer II',
                      'PL' => 'ILLUMINA',
                      'run' => 'WTSI:5143',
                      'LB' => 'WGS:WTSI:938',
                      'SM' => 'f393bb05-47c6-f3a5-e040-11ac0d48452c',
                      'type' => 'WGS',
                      'PI' => '480',
                      'md5' => 'PD4115b/5143_6.bam.md5',
                      'DT' => '2010-08-27T00:00:00.000+00:00'
                    };
  my $bam_ob_norm2 = {'exp' => 'WTSI:5143:938',
                      'info' => {
                                  'dcc_project_code' => 'BRCA-UK',
                                  'submitter_donor_id' => 'CGP_donor_1199137',
                                  'submitter_sample_id' => 'PD4115bx',
                                  'dcc_specimen_type' => 'Normal - blood derived',
                                  'use_cntl' => 'N/A',
                                  'total_lanes' => 21,
                                  'submitter_specimen_id' => 'CGP_specimen_1227970x',
                                  'ega_sample_accession' => 'EGAN00001000593x'
                                },
                      'CN' => 'WTSI',
                      'ID' => 'WTSI6066x',
                      'PU' => 'WTSI:5143_6x',
                      'file' => 'PD4115b/5143_6x.bam',
                      'PM' => 'Illumina Genome Analyzer II',
                      'PL' => 'ILLUMINA',
                      'run' => 'WTSI:5143x',
                      'LB' => 'WGS:WTSI:938x',
                      'SM' => 'f393bb05-47c6-f3a5-e040-11ac0d48452x',
                      'type' => 'WGS',
                      'PI' => '480',
                      'md5' => 'PD4115b/5143_6x.bam.md5',
                      'DT' => '2010-08-27T00:00:00.000+00:00'
                    };
  if($type eq 'tumour') {
    return $bam_ob_tum;
  }
  elsif($type eq 'normal') {
    return $bam_ob_norm;
  }
  elsif($type eq 'normal2') {
    return $bam_ob_norm2;
  }
  else {
    die "bam_ob data function doesn't know of a type '$type'\n";
  }
}
