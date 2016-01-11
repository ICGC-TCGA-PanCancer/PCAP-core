# this is a catch all to ensure all modules do compile
# added as lots of 'use' functionality is dynamic in pipeline
# and need to be sure that all modules compile.
# simple 'perl -c' is unlikely to work on head scripts any more.

use strict;
use Test::More;
use File::Which qw(which);
use List::Util qw(first);
use Const::Fast qw(const);
use Capture::Tiny qw(capture);
use Data::Dumper;
use version 0.77;

const my @REQUIRED_PROGRAMS => qw(bamcollate2 bammarkduplicates bamsort bwa);
const my $BIOBAMBAM2_VERSION => '2.0.25';
const my $BWA_VERSION => '0.7.12';

# can't put regex in const
my %EXPECTED_VERSION = (
                        'bamcollate2'       => {
                              'get'   => q{ -h},
                              'match' => qr/This is biobambam version ([[:digit:]\.]+)\./,
                              'version'       => version->parse($BIOBAMBAM2_VERSION)},
                        'bammarkduplicates' => {
                              'get'   => q{ -h},
                              'match' => qr/This is biobambam version ([[:digit:]\.]+)\./,
                              'version'       => version->parse($BIOBAMBAM2_VERSION)},
                        'bamsort'           => {
                              'get'   => q{ -h},
                              'match' => qr/This is biobambam version ([[:digit:]\.]+)\./,
                              'version'       => version->parse($BIOBAMBAM2_VERSION)},
                        'bwa'           => {
                              'get'   => q{},
                              'match' => qr/Version: ([[:digit:]\.]+[[:alpha:]]?)/, # we don't care about the revision number
                              'version'       => version->parse($BWA_VERSION)},
                        );

subtest 'External programs exist on PATH' => sub {
  for my $prog(@REQUIRED_PROGRAMS) {
    my $path = which($prog);
    isnt($path, q{}, "$prog found at $path");
  }
};

subtest 'External programs have expected version' => sub {
  for my $prog(@REQUIRED_PROGRAMS) {
    my $path = which($prog);
    my $details = $EXPECTED_VERSION{$prog};
    my $command = $path.$details->{'get'};
    my ($stdout, $stderr, $exit) = capture{ system($command); };
    my $reg = $details->{'match'};
    my ($version) = $stderr =~ /$reg/m;
    version->parse($version);

    ok(version->parse($version) >= $details->{'version'}, sprintf 'Expect minimum version of %s for %s', $details->{'version'}, $prog);
  }
};

done_testing();
