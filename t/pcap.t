use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Temp qw(tempdir);

const my $MODULE => 'PCAP';

my $obj;
subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};


ok(PCAP::license(), 'License text retrieved');

done_testing();

