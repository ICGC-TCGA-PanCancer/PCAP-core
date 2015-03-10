#!/usr/bin/perl

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014-2015 ICGC PanCancer Project
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not see:
#   http://www.gnu.org/licenses/gpl-2.0.html
##########LICENCE##########

use strict;
use warnings FATAL => 'all';
use English qw( -no_match_vars );
use Carp qw(croak);
use Getopt::Long;
use Pod::Usage;

use autodie;
use Capture::Tiny qw(:all);
use Config::IniFiles;
use Const::Fast qw(const);
use File::Copy qw(move);
use File::Fetch;
use File::Path qw(make_path);
use File::Which qw(which);
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use JSON qw(decode_json);
use List::Util qw(first any);


use Data::Dumper;

use PCAP;
use PCAP::Cli;

# perl-5.16.3 -I lib bin/gnos_pull.pl -o $SCRATCH112/gnos_pull -a CALLS -c ~/GitHub/PCAP-core/examples/cgp_gnos_pull.ini
# perl-5.16.3 -I lib bin/gnos_pull.pl -o $SCRATCH112/gnos_pull -a ALIGNMENTS -c ~/GitHub/PCAP-core/examples/cgp_gnos_pull.ini


const my @ANALYSIS_TYPES => (qw(ALIGNMENTS CALLS));
const my $DEFAULT_URL => 'http://pancancer.info/gnos_metadata/latest';
const my $GTDL_COMMAND => '%s -v -c %s -d %scghub/data/analysis/download/%s -p %s';

{
  my $options = option_builder();
  load_config($options);
  my $gz_manifest = manifest($options);
  my $to_process = load_data($options, $gz_manifest);
  pull_data($options, $to_process);
}

sub load_config {
  my $options = shift;
  my $cfg = Config::IniFiles->new( -file => $options->{'config'} );

  # load up filters, named by top level key in jsonl doc
  for my $param($cfg->Parameters('FILTERS')) {
    my @vals = $cfg->val('FILTERS', $param);
    next if(scalar @vals == 0);
    if(scalar @vals == 1) {
      my $val = shift @vals;
      next unless(defined $val && (length $val) > 0);
      @vals = split /,/, $val;
    }
    $options->{'FILTERS'}->{$param} = \@vals;
  }

  croak sprintf q{'KEY_FILE' Ssection is absent from %s}, $options->{'config'} unless($cfg->SectionExists('KEY_FILES'));
  for my $params($cfg->Parameters('KEY_FILES')) {
    $options->{'keys'}->{$params} = $cfg->val('KEY_FILES', $params);
  }

  if($cfg->exists('GENERAL', 'gtbin')) {
    $options->{'gtdownload'} = $cfg->val('GENERAL', 'gtbin').'/gtdownload';
  }
  else {
    $options->{'gtdownload'} = which('gtdownload');
  }

  if($cfg->exists('TRANSFER', 'order')) {
    $options->{'transfer'} = [$cfg->val('TRANSFER', 'order')]
  }

  return 1;
}

sub pull_data {
  my ($options, $to_process) = @_;

  my $outbase = $options->{'outdir'}.'/original';
  my $symbase = $options->{'outdir'}.'/donors';
  make_path($outbase) unless(-e $outbase);
  make_path($symbase) unless(-e $symbase);

  for my $donor(@{$to_process}) {
    my $donor_uniq = $donor->{'donor_unique_id'};
    my $donor_base = "$symbase/$donor_uniq";
    my $orig_base = "$outbase/$donor_uniq";
    make_path($orig_base) unless(-e $orig_base);
    pull_calls($options, $donor, $orig_base, $donor_base) if($options->{'analysis'} eq 'CALLS');
    pull_alignments($options, $donor, $orig_base, $donor_base) if($options->{'analysis'} eq 'ALIGNMENTS');
  }
}

sub pull_alignments {
  my ($options, $donor, $outbase, $donor_base) = @_;

  # for normal:
  pull_bam($options, $donor->{'normal_alignment_status'}, $outbase, $donor_base.'/normal');

  # for tumour
  for my $tumour_data(@{$donor->{'tumor_alignment_status'}}) {
    pull_bam($options, $tumour_data, $outbase, $donor_base.'/tumour');
  }
}

sub pull_bam {
  my ($options, $bam_data, $outbase, $sample_base) = @_;

  my $repo = select_repo($options, $bam_data->{'aligned_bam'}->{'gnos_repo'});
  my $gnos_id = $bam_data->{'aligned_bam'}->{'gnos_id'};

  my $download = sprintf $GTDL_COMMAND, $options->{'gtdownload'},
                                          $options->{'keys'}->{$repo},
                                          $repo,
                                          $gnos_id,
                                          $outbase;
  my $f_base = $outbase.'/'.$gnos_id;
  my $success = $f_base.'/.success.t';
  return if(-e $success);

  my $out_fh = IO::File->new("$f_base.out.log", "w+");
  my $err_fh = IO::File->new("$f_base.err.log", "w+");
  warn "Executing: $download\n";
  my $exit;
  capture { $exit = system($download); } stdout => $out_fh, stderr => $err_fh;
  $out_fh->close;
  $err_fh->close;
  die "A problem occured while executing: $download\n\tPlease check $f_base.err.log and proxy settings\n" if($exit);
  # as successful clean up logs
  unlink "$f_base.out.log";
  unlink "$f_base.err.log";

  make_path($sample_base) unless(-e $sample_base);
  my $aliq_id = $bam_data->{'aliquot_id'};
  my $orig_bam = $f_base.'/'.$bam_data->{'aligned_bam'}->{'bam_file_name'};
  my $sym_bam = $sample_base.'/'.$bam_data->{'aliquot_id'}.'.bam';
  symlink($orig_bam, $sym_bam) unless(-e $sym_bam);
  symlink($orig_bam.'.bai', $sym_bam.'.bai') unless(-e $sym_bam.'.bai');

  # touch a success file in the output loc
  open my $S, '>', $success; close $S;
}

sub pull_calls {
  my ($options, $donor, $outbase, $donor_base) = @_;
  # this should loop over the data found in 'is_sanger_vcf_in_jamboree' once modified to be array
  for my $caller(qw(sanger)) {
    # if we can get access to transfer metrics we want to select the fastest, or use list in config file

    my $repo = select_repo($options, $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_repo'});
    my $gnos_id = $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_id'};
    my $download = sprintf $GTDL_COMMAND, $options->{'gtdownload'},
                                          $options->{'keys'}->{$repo},
                                          $repo,
                                          $gnos_id,
                                          $outbase;
    my $f_base = $outbase.'/'.$gnos_id;
    my $success = $f_base.'/.success.t';
    next if(-e $success);

    my $out_fh = IO::File->new("$f_base.out.log", "w+");
    my $err_fh = IO::File->new("$f_base.err.log", "w+");
    warn "Executing: $download\n";
    my $exit;
    capture { $exit = system($download); } stdout => $out_fh, stderr => $err_fh;
    $out_fh->close;
    $err_fh->close;
    die "A problem occured while executing: $download\n\tPlease check $f_base.err.log and proxy settings\n" if($exit);
    # as successful clean up logs
    unlink "$f_base.out.log";
    unlink "$f_base.err.log";

    # now build symlinks
    make_path($donor_base) unless(-e $donor_base);
    symlink($f_base, "$donor_base/$caller") unless(-e "$donor_base/$caller");

    # touch a success file in the output loc
    open my $S, '>', $success; close $S;
  }
}

sub select_repo {
  my ($options, $available_at) = @_;
  return $available_at->[0] if(scalar @{$available_at} == 1);
  my $repo = $available_at->[0];
  for my $ordered_repo(@{$options->{'transfer'}}) {
    if(any {$_ eq $ordered_repo} @{$available_at}) {
      $repo = $ordered_repo;
      last;
    }
  }
  unless($repo eq $available_at->[0]) {
    warn "Changed repo from $available_at->[0] to $repo\n" if($options->{'debug'});
  }
  return $repo;
}

sub load_data {
  my ($options, $gz_manifest) = @_;
  my $total = 0;
  my @to_process;
  my @filter_keys = keys %{$options->{'FILTERS'}};
  my $z = IO::Uncompress::Gunzip->new($gz_manifest, MultiStream => 1) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  DONOR: while(my $jsonl = <$z>) {
    my $donor = decode_json($jsonl);

    $total++;

    # things to silently skip
    next if($donor->{'flags'}->{'is_test'});

    # user defined filters
    for my $key(@filter_keys) {
      unless(any{$_ eq $donor->{$key}} @{$options->{'FILTERS'}->{$key}}) {
        warn sprintf "Donor: %s not member of %s listed in %s\n", $donor->{'donor_unique_id'}, $key, $options->{'config'} if($options->{'debug'});
        next DONOR;
      }
    }

    # warn if debug on
    if($donor->{'flags'}->{'is_donor_blacklisted'}) {
      warn sprintf "Donor: %s blacklisted\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }

    if($donor->{'flags'}->{'is_manual_qc_failed'}) {
      warn sprintf "Donor: %s blacklisted\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }
    unless($donor->{'flags'}->{'is_normal_specimen_aligned'}) {
      warn sprintf "Donor: %s normal sample not aligned\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }
    unless($donor->{'flags'}->{'are_all_tumor_specimens_aligned'}) {
      warn sprintf "Donor: %s all tumour samples not aligned\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }

    if($options->{'analysis'} eq 'CALLS') {
      if(!$donor->{'flags'}->{'is_sanger_variant_calling_performed'} && $donor->{'flags'}->{'is_sanger_vcf_in_jamboree'}) {
        warn sprintf "Donor: %s is_sanger_variant_calling_performed=false BUT is_sanger_vcf_in_jamboree=true\n", $donor->{'donor_unique_id'} if($options->{'debug'});
        next;
      }

      unless($donor->{'flags'}->{'is_sanger_vcf_in_jamboree'}) {
        warn sprintf "Donor: %s results not in jamboree\n", $donor->{'donor_unique_id'} if($options->{'debug'});
        next;
      }
    }

    push @to_process, $donor;
  }
  close $z;
  my $retained = scalar @to_process;
  warn sprintf "Retained donors: %s\nRejected donors: %s\n", $retained, $total - $retained;
  return \@to_process;
}

sub manifest {
  my $options = shift;
  my $url = $options->{'url'} ? $options->{'url'} : $DEFAULT_URL;
  my $ff = File::Fetch->new(uri => $url);
  my $listing;
  my $where = $ff->fetch( to => \$listing );
  my ($file) = $listing =~ m/(donor_p_[[:digit:]]+[.]jsonl[.]gz)/xms;
  my $to_get = "$url/$file";
  make_path($options->{'outdir'}) unless(-e $options->{'outdir'});
  $ff = File::Fetch->new(uri => $to_get);
  $where = $ff->fetch( to => $options->{'outdir'});
  my $dest = $options->{'outdir'}.'/donor_p.jsonl.gz';
  move($where, $dest);
  return $dest;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'm|man' => \$opts{'m'},
		'u|url=s' => \$opts{'url'},
		'a|analysis=s' => \$opts{'analysis'},
		'c|config=s' => \$opts{'config'},
		'o|outdir=s' => \$opts{'outdir'},
		'd|debug' => \$opts{'debug'},
    'n|noproxy' => \$opts{'noproxy'},
    'v|version' => \$opts{'version'},
	) or pod2usage(2);

	pod2usage(-message => PCAP::license, -verbose => 0) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'m'});

  if(exists $opts{'version'} && defined $opts{'version'}) {
    print sprintf 'VERSION: %s\n', ();
    exit 0;
  }

  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
  PCAP::Cli::file_for_reading('config', $opts{'config'});
  pod2usage(-msg  => "\nERROR: Option '-a' is missing or invalid.\n", -verbose => 2,  -output => \*STDERR) unless(defined $opts{'analysis'} && first {$opts{'analysis'} eq $_} @ANALYSIS_TYPES);


	return \%opts;
}

__END__

=head1 NAME

gnos_pull.pl - retrieve/update analysis flow results on local systems.

=head1 SYNOPSIS

gnos_pull.pl [-h] -u http://XXX/gnos_metadata/latest/ -c pem_map.ini -o local_mirror/

  Required input:

    --analysis  (-a)  ALIGNMENTS or CALLS

    --outdir    (-o)  Where to save jsonl and resulting GNOS downloads

    --config    (-c)  Mapping of GNOS repos to permissions keys

  Other options:

    --url       (-u)  The base URL to retrieve jsonl file from
                        [http://pancancer.info/gnos_metadata/latest/]

    --debug     (-d)  prints extra debug information


    --help      (-h)  Brief documentation

    --man       (-m)  More verbose usage info

=head1 OPTIONS

=over 2

=item a|analysis

What type of data should be retrieved in this run.

ALIGNMENTS retrieves the mapped BAM for all samples linked to this donor.  Only retrieves data if all BAMs are ready.

CALLS retrieves the data for all variant calling workflows that have been run against the selected donors.
If data is available in subsequent executions for a new caller then this will be retrieved without affecting existing results.

=item o|outdir

The base output directory, this will need to be very large if pulling BAM files.

=item c|config

The config.ini file to be used.  There is a full example of this in the distribution.  Here is a brief example:

 [GENERAL]
 # used in preference to path, optional
 gtbin=/software/genetorrent/bin

 [FILTERS]
 # filters are applied as a union, so all must be fulfilled
 # key must match a top level key in jsonl file
 # csv:
 dcc_project_code=BRCA-UK,BRCA-US

 # or multiline:
 donor_unique_id=<<EOT
 BRCA-UK::CGP_donor_1167078
 BRCA-UK::CGP_donor_1187030
 EOT

 gnos_study=

 [KEY_FILES]
 https://gtrepo-ebi.annailabs.com/=/fullpath/XXX.pem
 https://gtrepo-etri.annailabs.com/=/fullpath/XXX.pem
 https://gtrepo-dkfz.annailabs.com/=/fullpath/XXX.pem

 [TRANSFER]
 # When data at multiple sites use this list to preferentially choose source
 # when not found in list just try the one found in the json
 order=<<EOT
 https://gtrepo-dkfz.annailabs.com/
 https://gtrepo-ebi.annailabs.com/
 https://gtrepo-etri.annailabs.com/
 EOT

=item u|url

Provide if you don't want to use the 'latest' jsonl file.

Recommended if you want to complete a first pass on a fixed file.

=item d|debug

Turns on debug information which can aid in diagnosing why a particular donor is not retrieved.

=back

=head1 METHODS

=head2 load_config

Loads the provided *.ini file provided to '-c'.  Embeds relevant data into options hash.

=head2 pull_data

Generic control method for retrieving data for a donor from GNOS, see pull_alignments and pull_calls.

=head2 pull_alignments

Abstracts calls to pull_bam to give a single method call which will retrieve all alignements associated with this donor.

=head2 pull_bam

Performs GNOS retrieval for a single BAM file (and index) to the 'original/<DONOR>/' output area.

Generates symlinks in the 'donors' output area to simplify use.

A call to the related pull_alignments method for a donor with 2 tumour samples will result in the
following directory structure.

 .../donors/DONOR/tumour/X.bam[.bai]
 .../donors/DONOR/tumour/Y.bam[.bai]
 .../donors/DONOR/normal/Z.bam[.bai]

=head2 pull_calls

Performs GNOS retrieval for all available variant calling.

Generates symlinks in the 'donors' output area to simplify use, which would look like:

  .../donors/DONOR/sanger/TUMOUR.*
  .../donors/DONOR/broad/TUMOUR.*
  ...

=head2 select_repo

If data is available from multiple GNOS repositories and TRANSFER.order is defined in the config.ini file
this method selects the repo highest in the list (if not found defaults to first in data structure).

NOTE - if up to date transfer metric reports become available these should be used in preference.

=head2 load_data

Parse the jsonl manifest, filter appropriately:

 Silently skip 'is_test' data.

If debug warn and skip when:

 Indicated by filters in *.ini
 Blacklisted
 Tumour or normal not aligned

If -a 'CALLS' skip when:

 No variant calling performed
 Not in Jamboree

=head2 manifest

Retrieve the donor_pXXX.jsonl.gz from pancancer.info or if provided alternate source.

=head2 option_builder

Parses the command line args and does basic checking of inputs.

=cut


