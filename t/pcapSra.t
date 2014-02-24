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
  my $tmp_pa = "$dir/cv_tables/TCGA/portionAnalyte.txt";
  unlink $tmp_pa;

  like(exception{PCAP::SRA::create_cv_lookups($dir)}
      , qr/ERROR: Unable to find controlled vocabulary file for /
      , 'Die when CV file not found');

  # re-create the missing file but empty
  `touch $tmp_pa`;

  like(exception{PCAP::SRA::create_cv_lookups($dir)}
      , qr/ERROR: Controlled vocabulary file for analyte_code is empty: /
      , 'Die when CV file is empty');
};

my $uuid = PCAP::SRA::uuid();
is($uuid, lc $uuid, 'Confirm uuid has been created in lower-case');

ok(PCAP::SRA::validate_seq_type('WGS'), 'Expected sequencing type');
like(exception{PCAP::SRA::validate_seq_type('XXX')}
  , qr/ERROR: '.*' is not a recognised sequencing\/library type/
  , 'Die when unexpected sequencing type');

done_testing();
