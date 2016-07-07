use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Temp qw(tempdir);

const my $MODULE => 'PCAP::Threaded';

my $obj;
subtest 'Initialisation checks' => sub {
  local $SIG{__WARN__}=sub{};
  use_ok($MODULE);
  like(exception{$MODULE->new('x')}, qr/Number of threads was NAN: /m, 'Requires integer thread count.');
  $obj = new_ok($MODULE);
  $obj = new_ok($MODULE => [1]);
};

subtest 'add_function checks' => sub {
  ok($obj->add_function('add_one', \&add_one), 'Added function add_one');
  like(exception{$obj->add_function('add_one', \&add_one)}
      , qr/Function add_one has already been defined./m
      , 'Fail when attempt to redefine function');
  like(exception{$obj->add_function('pass_arrref', [1])}
      , qr/Second argument to add_function should be a code reference, I got /m
      , 'Fail when not a coderef');
  like(exception{$obj->add_function('not_a_ref', &add_one)}
      , qr/Second argument to add_function should be a code reference, I got /m
      , 'Modify message when not a reference, where coderef expected');

};

subtest 'run checks' => sub {
  $obj = new_ok($MODULE => [1]);
  $obj->add_function('add_one', \&add_one);
  $obj->add_function('to_fail', \&to_fail);

  like(exception{$obj->run()}
      , qr/Iterations must be defined/
      , 'Fail when iterations not defined');
  like(exception{$obj->run('x')}
      , qr/Iterations must be a positive integer:/
      , 'Fail when iterations not a number');
  like(exception{$obj->run(0)}
      , qr/Iterations must be a positive integer:/
      , 'Fail when iterations == 0');
  like(exception{$obj->run(1)}
      , qr/Function_name must be defined/
      , 'Fail when function_name not defined');
  like(exception{$obj->run(1, 'non_existant')}
      , qr/Unable to find '.+', please check your declaration of add_function/
      , 'Fail when function not defined');
  ok($obj->run(1, 'add_one'), 'Success on clean sub');
  like(exception{$obj->run(1, 'to_fail')}
      , qr/Expected to fail/
      , 'Fail when function throws error');
## breaks Devel::Cover
#  $obj = new_ok($MODULE => [2]);
#  $obj->add_function('add_one', \&add_one);
#  ok($obj->run(2, 'add_one'), 'Success on multiple interations');
};

subtest 'thread object coniguration' => sub {
  is(&PCAP::Threaded::use_out_err, 1, 'Default value for out_err = 1');
  is(&PCAP::Threaded::disable_out_err, 0, 'Disabling out_err returns 0');
  is(&PCAP::Threaded::use_out_err, 0, 'Value following disable_out_err 0');
  is(&PCAP::Threaded::enable_out_err, 1, 'Enabling out_err returns 1');
  is(&PCAP::Threaded::use_out_err, 1, 'Value following enable_out_err 1');
  $obj = new_ok($MODULE => [1]);
  is($obj->thread_join_interval, 1, 'Default value for thread_join_interval = 1');
  is($obj->thread_join_interval(2), 2, 'Changing value for thread_join_interval = 2');
  is($obj->thread_join_interval, 2, 'Following change to 2, thread_join_interval = 2');
};

subtest 'completion utility checks' => sub {
  local $SIG{__WARN__}=sub{};
  my $dir = tempdir( CLEANUP => 1 );
  is(PCAP::Threaded::success_exists($dir, 1), 0, 'No success file');
  ok(PCAP::Threaded::touch_success($dir, 1), 'No success file');
  is(PCAP::Threaded::success_exists($dir, 1), 1, 'Success file present');
  ok(PCAP::Threaded::external_process_handler($dir, 'ls', 1), 'External process executes');
};

subtest 'thread divisor checks' => sub {
  like(exception{$obj->_suitable_threads('x')}
      , qr/Thread divisior must be a positive integer:/
      , 'Fail when divisor not a number');
  like(exception{$obj->_suitable_threads(0)}
      , qr/Thread divisior must be a positive integer:/
      , 'Fail when divisor == 0');
  is($obj->_suitable_threads(), 1, 'Return total_threads when no divisor');
  is($obj->_suitable_threads(2), 1, 'Return 1 when result is < 1');
};

done_testing();


sub add_one {
  return 1;
}

sub to_fail {
  die 'Expected to fail';
}
