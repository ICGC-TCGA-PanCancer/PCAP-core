use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use File::Path qw(make_path remove_tree);
use Try::Tiny qw(try catch finally);
use Fcntl qw( :mode );

use constant MODULE => 'PCAP::Cli';

use_ok(MODULE);

use FindBin qw($Bin);
my $test_data = "$Bin/../testData";


subtest 'file_for_reading' => sub {
  is(exception{ PCAP::Cli::file_for_reading('test', undef) }
      , qq{Option 'test' has not been defined.\n}
      , 'Fail when no filepath provided');
  like(exception{ PCAP::Cli::file_for_reading('test', File::Spec->catfile($test_data, 'no.file')) }
      , qr/Option 'test' \(.+\) should be an existing file./m
      , 'Fail when file found');
  like(exception{ PCAP::Cli::file_for_reading('test', File::Spec->catfile($test_data, 'empty.file')) }
      , qr/Option 'test' \(.+\) should be file with non-zero size./m
      , 'Fail when empty file found');
  like(exception{ PCAP::Cli::file_for_reading('test', $test_data) }
      , qr/Option 'test' \(.+\) is a directory./m
      , 'Fail when provided a directory');


  is(PCAP::Cli::file_for_reading('test', File::Spec->catfile($test_data, 'data.file')), File::Spec->catfile($test_data, 'data.file'), 'Pass when file with size found');
};

subtest 'out_dir_check' => sub {

  is(exception{ PCAP::Cli::out_dir_check('test', undef) }
      , qq{Option 'test' has not been defined.\n}
      , 'Fail when no path provided');
  like(exception{ PCAP::Cli::out_dir_check('test', File::Spec->catfile($test_data, 'data.file')) }
      , qr/Option 'test' points to an existing entity \(not a directory\):/m
      , 'Fail when pointed to file');

  my $tmp_dir = File::Spec->catdir($test_data, 'test_folder');
  my $non_write = File::Spec->catfile($test_data, 'nonWritableDir');
  # need to ensure folder is removed
  try {
    make_path($non_write);
    chmod S_IRUSR, $non_write;
    like(exception{ PCAP::Cli::out_dir_check('test', $non_write) }
        , qr/Option 'test' points to an existing WRITE PROTECTED directory:/m
        , 'Fail when pointed to non-writable area');
    is(PCAP::Cli::out_dir_check('test', $tmp_dir), $tmp_dir, 'Success when able to create directory');
    is(PCAP::Cli::out_dir_check('test', $tmp_dir), $tmp_dir, 'Success when directory exists and writable');
    ok((chmod 0400, $tmp_dir), 'make test folder readonly for next test');
    like(exception{ PCAP::Cli::out_dir_check('test', $tmp_dir) }
        , qr/Option '.+' points to an existing WRITE PROTECTED directory: /
        , 'Fail when provided an existing write protected dir.');
    my $unwriteable = File::Spec->catdir($tmp_dir, 'never_going_to_happen');
    like(exception{ PCAP::Cli::out_dir_check('test', $unwriteable) }
        , qr/Permission denied/
        , 'Fails when unable to create directory (parent protected)');

  } catch{ }
  finally {
    if(-e $tmp_dir) {
      chmod S_IRWXU, $tmp_dir; # if it don't work it's knackered anyway
      remove_tree($tmp_dir);
    }
    if(-e $non_write) {
      chmod S_IRWXU, $non_write;
      remove_tree($non_write);
    }
  };
};


done_testing();
