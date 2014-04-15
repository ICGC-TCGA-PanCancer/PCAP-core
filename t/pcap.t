use strict;
use Test::More;
use Const::Fast qw(const);

const my $MODULE => 'PCAP';

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};


ok(PCAP::license(), 'License text retrieved');

is(PCAP::upgrade_path(), 'biobambam,samtools,bwa', 'Default program install when no previous version');
is(PCAP::upgrade_path('9.9.9'), 'biobambam,samtools,bwa', 'Default program install when unknown version installed');

done_testing();

