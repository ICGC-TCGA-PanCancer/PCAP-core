use strict;
use Test::More;
use Test::Fatal;
use Const::Fast qw(const);
use File::Spec;
use FindBin qw($Bin);
use File::Temp qw(tempdir);

my $test_data = "$Bin/../testData";

const my $MODULE => 'PCAP::Bwa::Meta';
const my $REF_INIT => { 'in'    => 'somefile',
                        'temp' => 'somepath',};
const my $SET_RG_VAL => 5;
const my $RG_DEFAULT => qr/\@RG\tID:[a-z0-9]{8}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{12}\tCN:SANGER\tDS:short\tLB:SAMPLE_LIBRARY\tPI:500\tPL:HiSeq\tPU:1_1\tSM:SAMPLE_NAME/;
const my $RG_TAGS => {'SM' => 'wibble', 'LB' => 'wobble', };
const my $RG_STRING => qr/\@RG\\tID:[a-z0-9]{8}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{12}\\tCN:SANGER\\tDS:short\\tLB:wobble\\tPI:500\\tPL:HiSeq\\tPU:1_1\\tSM:wibble/;
const my $RG_PRINT => qr/\@RG\tID:[a-z0-9]{8}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{12}\tCN:SANGER\tDS:short\tLB:wobble\tPI:500\tPL:HiSeq\tPU:1_1\tSM:wibble/;
const my @VALID_FASTQ_EXT => qw(fastq fq fastq.gz fq.gz);


#my $fastq_paired_one = File::Spec->catfile($test_data, 'test_1.fastq');

my $bail_out = 0;

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
  my $obj = new_ok($MODULE => [$REF_INIT]);
};

subtest 'Non-object funcions' => sub {

  subtest 'FASTQ: extensions' => sub {
    for my $test_ext(@VALID_FASTQ_EXT) {
      my $test_file = "file.$test_ext";
      my ($file, $ext) = PCAP::Bwa::Meta::is_fastq_ext($test_file);
      isnt($file, $test_file, "$test_ext file should have extension removed") or $bail_out++;
      is($ext, $test_ext, "Expect extension $test_ext") or $bail_out++;
    }
  };

  BAIL_OUT('Subsequent fastq tests will fail') if($bail_out);

  subtest 'bam extension' => sub {
    my $bamfile = 'file.bam';
    my ($file, $ext) = PCAP::Bwa::Meta::is_fastq_ext($bamfile);
    is($file, $bamfile, "bam file should not have extension removed");
    ok(!defined $ext, "Extension should not be set");
  };

  subtest 'FASTQ: interleaved or paired' => sub {
    my ($base, $end);
    ($base, $end) = PCAP::Bwa::Meta::parse_fastq_filename('fastq_1');
    is($base, 'fastq', "fastq_1 file should have '_1' removed");
    is($end, '1', "Read end should be 1");
    ($base, $end) = PCAP::Bwa::Meta::parse_fastq_filename('fastq_2');
    is($base, 'fastq', "fastq_2 file should have '_2' removed");
    is($end, '2', "Read end should be 2");
    ($base, $end) = PCAP::Bwa::Meta::parse_fastq_filename('fastq_3');
    isnt($base, 'fastq', "fastq_3 file should be unchanged, assume interleaved");
    ok(!defined $end, "_3, so not defined");
    ($base, $end) = PCAP::Bwa::Meta::parse_fastq_filename('fastq');
    is($base, 'fastq', "fastq file should be unchanged, assume interleaved");
    ok(!defined $end, "no end trailer, so not defined");
  };

  subtest 'set_rg_index inputs' => sub {
    is(PCAP::Bwa::Meta::set_rg_index(5), -1, 'Set rg_index to positive integer passes');
    is(PCAP::Bwa::Meta::set_rg_index(15), -1, 'Set rg_index to pos-larger integer passes');

    like( exception{ PCAP::Bwa::Meta::set_rg_index() }
        , qr/set_rg_index requires a value/
        , 'Set rg_index to float fails');
    like( exception{ PCAP::Bwa::Meta::set_rg_index(1.5) }
        , qr/Value must be a positive integer/
        , 'Set rg_index to float fails');
    like( exception{ PCAP::Bwa::Meta::set_rg_index('wibble') }
        , qr/Value must be a positive integer/
        , 'Set rg_index to string fails');
    like( exception{ PCAP::Bwa::Meta::set_rg_index(0) }
        , qr/Value must be a positive integer/
        , 'Set rg_index to 0 fails');
    like( exception{ PCAP::Bwa::Meta::set_rg_index(-1) }
        , qr/Value must be a positive integer/
        , 'Set rg_index to -1 fails');
  };
};

subtest 'Objects from file list' => sub {
  my $tmp = tempdir( CLEANUP => 1 );
  is(&PCAP::Bwa::Meta::reset_rg_index, -1, q{Reset of rg_index}); # the rest will fail if this hasn't worked
  like( exception {PCAP::Bwa::Meta::files_to_meta() }
      , qr/Requires tmpdir and array-ref of files/
      , 'Check inputs present');
  like( exception {PCAP::Bwa::Meta::files_to_meta(undef, []) }
      , qr/Requires tmpdir and array-ref of files/
      , 'Check inputs present, missing $tmp');
  like( exception {PCAP::Bwa::Meta::files_to_meta($tmp) }
      , qr/Requires tmpdir and array-ref of files/
      , 'Check inputs present, missing $files');
  like( exception {PCAP::Bwa::Meta::files_to_meta('wibble', []) }
      , qr/Directory must exist: /
      , 'Check temp directory is present');
  like( exception {PCAP::Bwa::Meta::files_to_meta($tmp, {}) }
      , qr/\$files must be an array-ref/
      , 'Check $files is array-ref');
  like( exception {PCAP::Bwa::Meta::files_to_meta($tmp, []) }
      , qr/Some files must be provided/
      , 'Check $files is not empty array-ref');
  is(PCAP::Bwa::Meta::files_to_meta($tmp, [File::Spec->catfile($test_data, 'not_really_a.bam')])->[0]->{'fastq'}
    , undef, 'Bam not seen as fastq');
  is(PCAP::Bwa::Meta::files_to_meta($tmp, [File::Spec->catfile($test_data, '1.fq')], 'sample')->[0]->{'fastq'}
    , 'fq', 'Interleaved fq identified as fastq');
  is(PCAP::Bwa::Meta::files_to_meta($tmp, [File::Spec->catfile($test_data, '1.fq')], 'sample')->[0]->{'paired_fq'}
    , undef, 'Interleaved fq identified as not paired');
  is(PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, '1_1.fq')
                                          , File::Spec->catfile($test_data, '1_2.fq')], 'sample')->[0]->{'paired_fq'}
    , 1, 'Paired fq identified as paired');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp
                                                  , [ File::Spec->catfile($test_data, '2_1.fq')
                                                    , File::Spec->catfile($test_data, '2_2.fq')]
                                                  , 'sample') }
      , qr/Unable to find file for read 2, for /
      , 'Fail when file for second end missing');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp
                                                      , [ File::Spec->catfile($test_data, '3_1.fq')
                                                        , File::Spec->catfile($test_data, '3_2.fq')]
                                                      , 'sample')}
      , qr/Unable to find file for read 1, for /
      , 'Fail when file for first end missing');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp
                                                      , [ File::Spec->catfile($test_data, 'empty_r1_1.fq')
                                                        , File::Spec->catfile($test_data, 'empty_r1_2.fq')]
                                                      , 'sample') }
      , qr/File for read 1 is empty: /
      , 'Fail when file for first end is empty');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp
                                                      , [ File::Spec->catfile($test_data, 'empty_r2_1.fq')
                                                        , File::Spec->catfile($test_data, 'empty_r2_2.fq')]
                                                      , 'sample') }
      , qr/File for read 2 is empty: /
      , 'Fail when file for second end is empty');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, 'missing.fq')]
                                                      , 'sample') }
      , qr/File does not exist: /
      , 'Fail when interleaved fastq not found');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, 'empty.fq')]
                                                      , 'sample') }
      , qr/File is empty: /
      , 'Fail when interleaved fastq empty');
  is( exception{ PCAP::Bwa::Meta::files_to_meta($tmp
                                                    , [ File::Spec->catfile($test_data, '1_1.fq')
                                                      , File::Spec->catfile($test_data, '1_2.fq')
                                                      , File::Spec->catfile($test_data, '1.fq')]
                                                      , 'sample') }
      , qq{ERROR: BAM, paired FASTQ and interleaved FASTQ file types cannot be mixed, please choose one type\n}
      , 'Fail when inputs are mixed file types');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, 'wibble.wobble')]
                                                      , 'sample') }
      , qr/.+ is not an expected input file type.\n/m
      , 'Fail for unexpected filetype');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, 'missing.bam')]
                                                      , 'sample') }
      , qr/File does not exist: /
      , 'Fail when bam not found');
  like( exception{ PCAP::Bwa::Meta::files_to_meta($tmp, [ File::Spec->catfile($test_data, 'empty.bam')]
                                                      , 'sample') }
      , qr/File is empty: /
      , 'Fail when bam is empty');
};

my $meta; # for reuse
subtest 'Accessors' => sub {
  my $tmp = tempdir( CLEANUP => 1 );
  my $tstub = File::Spec->catfile($tmp, '1');
  is(&PCAP::Bwa::Meta::reset_rg_index, -1, q{Reset of rg_index}); # the rest will fail if this hasn't worked
  my $meta_set = PCAP::Bwa::Meta::files_to_meta($tmp, [File::Spec->catfile($test_data, '1_1.fq')], 'sample');
  $meta = $meta_set->[0];
  is($meta->in, File::Spec->catfile($test_data, '1'), 'Expected in');
  like( exception { $meta->in(1) }
      , qr/'in' can only be set via new\(\)/
      , 'Fail to directly set in');

  is($meta->fastq, 'fq', 'Expected fastq');
  like( exception { $meta->fastq(1) }
      , qr/'fastq' can only be set via new\(\)/
      , 'Fail to directly set fastq');

  is($meta->paired_fq, 1, 'Expected paired_fastq');
  like( exception { $meta->paired_fq(1) }
      , qr/'paired_fq' can only be set via new\(\)/
      , 'Fail to directly set paired_fq');

  is($meta->tstub, $tstub, 'Generate expected tstub');
  ok(exists $meta->{'tstub'}, 'Underlying tstub hash val has been created');
  is($meta->{'tstub'}, $tstub, 'Underlying tstub hash val expected');
  is($meta->tstub, $tstub, 'Recover tstub from stored val');
  like( exception { $meta->tstub(1) }
      , qr/'tstub' is autopopulated/
      , 'Fail to directly set tstub');

  like( exception { $meta->rg(1) }
      , qr/'rg' is autopopulated/
      , 'Fail to directly set rg');
};

subtest 'rg_header checks' => sub {
  $meta = new_ok($MODULE => [{ 'in'    => File::Spec->catfile($test_data, 'header.bam'),
                                          'temp' => 'somepath',}]);

  like($meta->rg_header(qq{\t}), $RG_DEFAULT, 'RG default header constructed correctly');

  $meta = new_ok($MODULE => [{ 'in'    => File::Spec->catfile($test_data, 'header.bam'),
                                          'temp' => 'somepath',}]);

  like($meta->rg_header(q{\t}, $RG_TAGS), $RG_STRING, 'RG header constructed correctly');

  like( exception { $meta->rg_header(q{\t}, $RG_TAGS) }
      , qr/'rg_header' has already been set/
      , 'Fail to set rg_header a second time');

  like($meta->rg_header(q{\t}), $RG_STRING, 'RG header retrieved for arg pass');
  like($meta->rg_header(qq{\t}), $RG_PRINT, 'RG header retrieved for print');

  # clear header for further tests
  my $tmp = tempdir( CLEANUP => 1 );
  my $meta_set = PCAP::Bwa::Meta::files_to_meta($tmp, [File::Spec->catfile($test_data, '1_1.fq')]);
  $meta = $meta_set->[0];
  like( exception{ $meta->rg_header('.') }
      , qr/'.+' is manditory for RG header/
      , 'Fail if manditory elements of RG are missing');
};

subtest 'init hash options' => sub {
  like( exception{ $meta->_init({'rg', => 6}) }
      , qr/'rg' is auto-populated, to initialise a start value see PCAP::Bwa::Meta::set_rg_index/
      , 'Throw error on attempt to set rg in init hash');
  like( exception{ $meta->_init({'wibble', => 6}) }
      , qr/'wibble' is not a valid parameter for object initialisation/
      , 'Throw error on attempt to set wibble in init hash');
  like( exception{ $meta->_init({'fastq', => [1]}) }
      , qr/'fastq' is not a scalar, only simple values are expected/
      , 'Throw error on attempt to set a reference as a value');
  like( exception{ $meta->_init({ }) }
      , qr/'.+' must be exist/
      , q{Throw error when required keys don't exist});
  like( exception{ $meta->_init({ 'in' => undef }) }
      , qr/'.+' must be defined/
      , q{Throw error when required keys don't have a defined value});
  like( exception{ $meta->_init({ 'in' => q{} }) }
      , qr/'.+' must have value with non-0 length/
      , q{Throw error when required keys don't have a value with non-0 length});
};

done_testing();
