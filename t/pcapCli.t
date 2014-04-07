use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use File::Path qw(make_path remove_tree);
use Try::Tiny qw(try catch finally);
use Fcntl qw( :mode );
use File::Temp;
use English qw( -no_match_vars );

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

  my $tmp_root = File::Temp->newdir('PCAPtests_XXXX');
  my $tmp_dir = File::Spec->catdir($tmp_root, 'test_folder');
  note("Temporary folder: $tmp_dir");
  # need to ensure folder is removed
  try {
    is(PCAP::Cli::out_dir_check('test', $tmp_dir), $tmp_dir, 'Success when able to create directory');
    is(PCAP::Cli::out_dir_check('test', $tmp_dir), $tmp_dir, 'Success when directory exists and writable');

    SKIP: {
      skip q{Running as superuser, certain file tests always true.}, 3 if($EFFECTIVE_USER_ID == 0);
      ok((chmod S_IRUSR, $tmp_dir), 'make test folder readonly for next test');
      like(exception{ PCAP::Cli::out_dir_check('test', $tmp_dir) }
          , qr/Option '.+' points to an existing WRITE PROTECTED directory: /
          , 'Fail when provided an existing write protected dir.');
      my $unwriteable = File::Spec->catdir($tmp_dir, 'never_going_to_happen');
      like(exception{ PCAP::Cli::out_dir_check('test', $unwriteable) }
          , qr/Permission denied/
          , 'Fails when unable to create directory (parent protected)');
    }
    SKIP: {
      skip q{Running as normal user, can't do superuser check.}, 1, if($EFFECTIVE_USER_ID != 0);
      like(exception{ PCAP::Cli::out_dir_check('test', $tmp_dir, 1) }
          , qr/EXIT: Please run as non-root user/
          , 'Fail when executed as super-user');
    }
  } catch{ }
  finally {
    if(-e $tmp_dir) {
      chmod S_IRWXU, $tmp_dir; # if it don't work it's knackered anyway
    }
  };
};


done_testing();
