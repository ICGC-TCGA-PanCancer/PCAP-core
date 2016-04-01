use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use Try::Tiny qw(try catch finally);
use English qw( -no_match_vars );
use Const::Fast;
use Capture::Tiny ':all';

use PCAP::Bam::Coverage;

use FindBin qw($Bin);
my $test_data = "$Bin/data";
my $script = "$Bin/../bin/xam_coverage_bins.pl";
my $input_bam = $test_data."/coverage.bam";
my $bed = $test_data."/coverage_exons.bed";
my $gff = $test_data."/coverage_exons.gff3";
my $test_output = $test_data."/test_cvg.txt";

const my $PERL => $^X;
const my $COMMAND_STDOUT => "$PERL %s -f %s -r %s -t %s";
const my $COMMAND_FILE => "$PERL %s -f %s -r %s -t %s -o %s";
const my $BED_TYPE => 'bed';
const my $GFF_TYPE => 'gff3';

const my $EXP_RESULT => '0:0,1-10:1.0000,11-20:1.0000,21-30:0.0000,31-40:0.0000,41-50:0.0000,51-100:0.0000,101-200:0.0000,201-500:0.0000,501+:0.0000';

subtest 'Object access BED' => sub {
  my $cov = new_ok('PCAP::Bam::Coverage', [{xam => $input_bam,
                                            type => $BED_TYPE,
                                            target => $bed}]);
  is($cov->build_depth, $EXP_RESULT, "Expected result BED");
};

subtest 'Object access BED' => sub {
  my $cov = new_ok('PCAP::Bam::Coverage', [{xam => $input_bam,
                                            type => $GFF_TYPE,
                                            target => $gff}]);
  is($cov->build_depth, $EXP_RESULT, "Expected result GFF3");
};


subtest 'bed no output file' => sub {
  my $cmd = sprintf($COMMAND_STDOUT,$script,$input_bam,$bed,$BED_TYPE);
  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };
  is($exit, 0,"Zero exit code for script");
  chomp($stdout);
  is($stdout, $EXP_RESULT,"Expected result to stdout");
};

subtest 'bed with output file' => sub {
  my $cmd = sprintf($COMMAND_FILE,$script,$input_bam,$bed,$BED_TYPE,$test_output);
  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };
  is($exit, 0,"Zero exit code for script");
  my $FH;
  ok(open($FH, '<', $test_output),"Read from test output file");
  my $line = <$FH>;
  ok(close($FH));
  chomp($line);
  is($line, $EXP_RESULT,"Expected result to stdout");
  ok(unlink($test_output),"Delete test output");
};

subtest 'gff3 no output file' => sub {
  my $cmd = sprintf($COMMAND_STDOUT,$script,$input_bam,$gff,$GFF_TYPE);
  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };
  is($exit, 0,"Zero exit code for script");
  chomp($stdout);
  is($stdout, $EXP_RESULT,"Expected result to stdout");
};

subtest 'gff3 with output file' => sub {
  my $cmd = sprintf($COMMAND_FILE,$script,$input_bam,$gff,$GFF_TYPE,$test_output);
  my ($stdout, $stderr, $exit) = capture {
    system( $cmd );
  };
  is($exit, 0,"Zero exit code for script");
  my $FH;
  ok(open($FH, '<', $test_output),"Read from test output file");
  my $line = <$FH>;
  ok(close($FH));
  chomp($line);
  is($line, $EXP_RESULT,"Expected result to stdout");
  ok(unlink($test_output),"Delete test output");
};





done_testing();
