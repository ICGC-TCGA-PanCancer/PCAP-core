use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);

const my $MODULE => 'PCAP::Bwa';

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'Non object checks' => sub {
  ok(PCAP::Bwa::bwa_version(), 'Version returned for BWA');
};

done_testing();
