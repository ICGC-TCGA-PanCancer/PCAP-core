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
use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw(croak);
use Getopt::Long;
use Pod::Usage;

use autodie;
use Capture::Tiny qw(:all);
use Config::IniFiles;
use Const::Fast qw(const);
use File::Copy qw(move);
use File::Fetch;
use File::Path qw(make_path remove_tree);
use File::Spec;
use IO::File;
use JSON qw(decode_json);
use List::Util qw(first any);
use Proc::PID::File;
use Data::Dumper;

use PCAP;
use PCAP::Cli;

const my @ANALYSIS_TYPES => (qw(ALIGNMENTS CALLS));
const my @AVAILABLE_COMPOSITE_FILTERS => (qw(not_sanger_workflow caller max_dataset_GB multi_tumour sanger_version broad_version dkfz_embl_version jamboree_approved manual_donor_blacklist));
const my $DEFAULT_URL => 'http://pancancer.info/gnos_metadata/latest';
const my $GTDL_COMMAND => '%s%s --max-children 3 --rate-limit 200 -vv -c %s -d %scghub/data/analysis/download/%s -p %s';

{
  my $options = option_builder();

  my $proc = Proc::PID::File->new(-verify => 1);
  $proc->file('dir' => $options->{'outdir'}, 'name' => 'gnos_pull.pl.'.$options->{'analysis'});
  if($proc->alive) {
    warn "Already running against output location: '$options->{outdir}' ... exiting!\n" if($options->{'debug'});
    exit 0;
  }
  $proc->touch;

  load_config($options);
  my $gz_manifest = manifest($options);
  my $to_process = load_data($options, $gz_manifest);
  exit 0 if(defined $options->{'info'});
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

  $options->{'COMPOSITE_FILTERS'}->{'jamboree_approved'} = 0;
  for my $comp_fil($cfg->Parameters('COMPOSITE_FILTERS')) {
    my $val = $cfg->val('COMPOSITE_FILTERS', $comp_fil);
    next unless($val);
    die sprintf "COMPOSITE_FILTERS.%s is not a supported filter type\n", $comp_fil unless(any{$comp_fil eq $_}@AVAILABLE_COMPOSITE_FILTERS);
    $options->{'COMPOSITE_FILTERS'}->{$comp_fil} = $val;
  }

  if(exists $options->{'COMPOSITE_FILTERS'}->{'manual_donor_blacklist'}) {
    my %bl_donors = map { $_ => 1 } split /\n/, $options->{'COMPOSITE_FILTERS'}->{'manual_donor_blacklist'};
    $options->{'COMPOSITE_FILTERS'}->{'manual_donor_blacklist'} = \%bl_donors;
  }
  if(exists $options->{'COMPOSITE_FILTERS'}->{'not_sanger_workflow'}) {
    my %bl_workflows = map { $_ => 1 } split /[\n,]/, $options->{'COMPOSITE_FILTERS'}->{'not_sanger_workflow'};
    $options->{'COMPOSITE_FILTERS'}->{'not_sanger_workflow'} = [keys \%bl_workflows];
  }

  croak sprintf q{'KEY_FILE' Ssection is absent from %s}, $options->{'config'} unless($cfg->SectionExists('KEY_FILES'));
  for my $params($cfg->Parameters('KEY_FILES')) {
    $options->{'keys'}->{$params} = $cfg->val('KEY_FILES', $params);
  }

  if($cfg->exists('GENERAL', 'gtbin')) {
    my $gtbin = $cfg->val('GENERAL', 'gtbin');
    $options->{'gtdownload'} = $gtbin.'/gtdownload';
    my $gtshare = $gtbin;
    $gtshare =~ s|bin/?$|share/GeneTorrent|;
    $options->{'gtdownload'} .= ' -R '.$gtshare;
  }
  else {
    $options->{'gtdownload'} = _which('gtdownload');
  }
  if($cfg->exists('GENERAL', 'gt_timeout_min')) {
    $options->{'gt_timeout'} = ' -k '.$cfg->val('GENERAL', 'gt_timeout_min');
  }
  else {
    $options->{'gt_timeout'} = q{};
  }


  if($cfg->exists('TRANSFER', 'order')) {
    $options->{'transfer'} = [$cfg->val('TRANSFER', 'order')]
  }

  return 1;
}

use threads;

sub pull_data {
  my ($options, $to_process) = @_;

  my $outbase = $options->{'outdir'}.'/original';
  my $symbase = $options->{'outdir'}.'/donors';
  make_path($outbase) unless(-e $outbase);
  make_path($symbase) unless(-e $symbase);

  my $code_ref;
  my $check_ref;
  if($options->{'analysis'} eq 'CALLS') {
    $code_ref = \&pull_calls;
  }
  elsif($options->{'analysis'} eq 'ALIGNMENTS') {
    $check_ref = \&check_alignments;
    $code_ref = \&pull_alignments;
  }

  my $thread_count = $options->{'threads'};


  while(@{$to_process} > 0) {
    my $donor = shift @{$to_process};
    my $donor_uniq = $donor->{'donor_unique_id'};
    # is PROJECT::donor, so convert to folder
    $donor_uniq =~ s|[:]{2}|/|;
    my $donor_base = "$symbase/$donor_uniq";
    my $orig_base = "$outbase/$donor_uniq";
    make_path($orig_base) unless(-e $orig_base);

    if(defined $check_ref) {
      next if($check_ref->($options, $donor, $orig_base, $donor_base) == 0);
    }

    warn "Submitting: $donor->{donor_unique_id}\n";
    if($thread_count > 1) {
      if(threads->list(threads::all) < $thread_count) {
        threads->create($code_ref, $options, $donor, $orig_base, $donor_base);
        # don't sleep if not full yet
        next if(threads->list(threads::all) < $thread_count);
      }
      sleep 1 while(threads->list(threads::joinable) == 0);
      for my $thr(threads->list(threads::joinable)) {
        $thr->join;
        if(my $err = $thr->error) { die "Thread error: $err\n"; }
      }
    }
    else {
      $code_ref->($options, $donor, $orig_base, $donor_base);
    }
  }
  if($thread_count > 1) {
    # last gasp for any remaining threads
    sleep 1 while(threads->list(threads::running) > 0);
    for my $thr(threads->list(threads::joinable)) {
      $thr->join;
      if(my $err = $thr->error) { die "Thread error: $err\n"; }
    }
  }
}

sub check_or_create_symlink {
  my ($source, $target) = @_;
  unless(-e $source) {
    die "FATAL: Source file doesn't exist: $source\n";
  }
  if(-l $target) {
    my $existing_target = readlink $target;
    if($source ne $existing_target) {
      unlink $target;
      symlink($source, $target);
    }
  }
  else {
    symlink($source, $target);
  }
  return 1;
}

sub check_alignments {
  my ($options, $donor, $outbase, $donor_base) = @_;
  warn "Checking $donor->{donor_unique_id}\n";
  my $to_do = 0;
  # for normal:
  $to_do += check_bam($options, $donor->{'donor_unique_id'}, $donor->{'normal_alignment_status'}, $outbase, $donor_base, 'normal');

  # for tumour
  for my $tumour_data(@{$donor->{'tumor_alignment_status'}}) {
    $to_do += check_bam($options, $donor->{'donor_unique_id'}, $tumour_data, $outbase, $donor_base, 'tumour');
  }
  return $to_do;
}

sub check_bam {
  my ($options, $donor_id, $bam_data, $outbase, $donor_base, $type) = @_;
  my $repo = select_repo($options, $bam_data->{'aligned_bam'}->{'gnos_repo'});
  unless(exists $options->{'keys'}->{$repo}) {
    warn sprintf "Skipping %s BAM for Donor %s - No permission key for repo %s", $type, $donor_id, $repo;
    return 0;
  }

  my $gnos_id = $bam_data->{'aligned_bam'}->{'gnos_id'};

  my $f_base = $outbase.'/'.$gnos_id;
  my $success = $f_base.'/.success.t';

  my $sample_base = "$donor_base/$type";

  make_path($sample_base) unless(-e $sample_base);
  my $aliq_id = $bam_data->{'aliquot_id'};
  my $orig_bam = $f_base.'/'.$bam_data->{'aligned_bam'}->{'bam_file_name'};
  my $sym_bam = $sample_base.'/'.$bam_data->{'aliquot_id'}.'.bam';

  if(-e $success) {
    check_or_create_symlink($orig_bam, $sym_bam);
    check_or_create_symlink($orig_bam.'.bai', $sym_bam.'.bai');
    create_bas($repo, $gnos_id, $sym_bam);
    return 0;
  }

  return 0 if($options->{'symlinks'});
  return 1;
}

sub pull_alignments {
  my ($options, $donor, $outbase, $donor_base) = @_;

  # for normal:
  pull_bam($options, $donor->{'donor_unique_id'}, $donor->{'normal_alignment_status'}, $outbase, $donor_base, 'normal');

  # for tumour
  for my $tumour_data(@{$donor->{'tumor_alignment_status'}}) {
    pull_bam($options, $donor->{'donor_unique_id'}, $tumour_data, $outbase, $donor_base, 'tumour');
  }
}

sub pull_bam {
  my ($options, $donor_id, $bam_data, $outbase, $donor_base, $type) = @_;

  my $repo = select_repo($options, $bam_data->{'aligned_bam'}->{'gnos_repo'});
  unless(exists $options->{'keys'}->{$repo}) {
    warn sprintf "Skipping %s BAM for Donor %s - No permission key for repo %s", $type, $donor_id, $repo;
    return 0;
  }
  my $gnos_id = $bam_data->{'aligned_bam'}->{'gnos_id'};

  my $f_base = $outbase.'/'.$gnos_id;
  my $success = $f_base.'/.success.t';

  my $sample_base = "$donor_base/$type";

  make_path($sample_base) unless(-e $sample_base);
  my $aliq_id = $bam_data->{'aliquot_id'};
  my $orig_bam = $f_base.'/'.$bam_data->{'aligned_bam'}->{'bam_file_name'};
  my $sym_bam = $sample_base.'/'.$bam_data->{'aliquot_id'}.'.bam';

  if(-e $success) {
    check_or_create_symlink($orig_bam, $sym_bam);
    check_or_create_symlink($orig_bam.'.bai', $sym_bam.'.bai');
    create_bas($repo, $gnos_id, $sym_bam);
    return;
  }
  return if($options->{'symlinks'});

  my $download = sprintf $GTDL_COMMAND, $options->{'gtdownload'},
                                        $options->{'gt_timeout'},
                                        $options->{'keys'}->{$repo},
                                        $repo,
                                        $gnos_id,
                                        $outbase;

  my $out_file = "$f_base.out.log";
  my $err_file = "$f_base.err.log";
  $download .= " > $f_base.out.log";
  $download = "($download) >& $err_file";
  warn "Executing: $download\n";

  my $exit = system($download);

  return 0 if($exit && download_abandoned($err_file));

  if($exit) {
    warn "ERROR: A problem occured while executing: $download\n\tPlease check $err_file and proxy settings\n";
    return 0;
  }
  # as successful clean up logs
  unlink $out_file;
  unlink $err_file;

  check_or_create_symlink($orig_bam, $sym_bam);
  check_or_create_symlink($orig_bam.'.bai', $sym_bam.'.bai');
  my $bas_valid = create_bas($repo, $gnos_id, $sym_bam);

  if($bas_valid == 1) {
    # touch a success file in the output loc
    open my $S, '>', $success; close $S;
  }
  return 1;
}

sub create_bas {
  my ($repo, $gnos_id, $sym_bam) = @_;
  my $bas_file = $sym_bam.'.bas';
  unlink $bas_file if(-l $bas_file);
  return 1 if(-s $bas_file);
  # pull the bas file, done here to handle back fill of this data
  my $get_bas = sprintf '%s %s/xml_to_bas.pl -d %scghub/metadata/analysisFull/%s -o %s -b %s',
                      $^X,
                      $Bin,
                      $repo,
                      $gnos_id,
                      $bas_file,
                      $sym_bam;
  my $success = 0;
  my $bas_valid = 1;
  for(0..5) {
    warn "Executing: $get_bas\n";
    my ($stdout, $stderr, $exit_c) = capture { system($get_bas); };
    if($exit_c) {
      if($stderr =~ m/Unable to recover read_group_id clash using PU field/ms) {
        warn "Aborting due to bad data: $get_bas\n\nSTDOUT:\n$stdout\n\nSTDERR:\n$stderr\n";
        $bas_valid = 0;
        last;
      }
      warn "A problem occured while executing: $get_bas\n\nSTDOUT:\n$stdout\n\nSTDERR:\n$stderr\n...Retry\n";
      sleep 30;
    }
    else {
      $success++;
      last;
    }
  }
  if($bas_valid) {
    die "Failed after multiple attempts, aborting bas generation using $get_bas" unless($success);
  }
  return $bas_valid;
}

sub pull_calls {
  my ($options, $donor, $outbase, $donor_base) = @_;
  my @callers = @{$donor->{'flags'}->{'variant_calling_performed'}};
  # this should loop over the data found in 'is_sanger_vcf_in_jamboree' once modified to be array
  for my $caller(@callers) {
    next if($caller eq 'embl' || $caller eq 'dkfz');

    next if(exists  $options->{'COMPOSITE_FILTERS'}->{'caller'} && index(','.$options->{'COMPOSITE_FILTERS'}->{'caller'}.',', ','.$caller.',') == -1);
    next unless($donor->{'flags'}->{'is_'.$caller.'_variant_calling_performed'});

    next if($options->{'COMPOSITE_FILTERS'}->{'jamboree_approved'} == 1 && !$donor->{'flags'}->{'is_'.$caller.'_vcf_in_jamboree'});

    next if(exists $options->{'COMPOSITE_FILTERS'}->{$caller.'_version'} && $options->{'COMPOSITE_FILTERS'}->{$caller.'_version'} ne $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'workflow_details'}->{'variant_workflow_version'});

    if(exists $options->{'COMPOSITE_FILTERS'}->{not_sanger_workflow} && first { $donor->{'variant_calling_results'}->{'sanger_variant_calling'}->{'workflow_details'}->{'variant_workflow_name'} eq $_ } @{$options->{'COMPOSITE_FILTERS'}->{not_sanger_workflow}}) {
      warn "Skipping version $donor->{donor_unique_id} as version == $options->{COMPOSITE_FILTERS}->{not_sanger_workflow}\n" if($options->{debug});
      next;
    }

    # if we can get access to transfer metrics we want to select the fastest, or use list in config file
    my $repo = select_repo($options, $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_repo'});
    unless(exists $options->{'keys'}->{$repo}) {
      warn sprintf "Skipping caller %s for Donor %s - No permission key for repo %s", $caller, $donor->{'donor_unique_id'}, $repo;
      next;
    }


    next unless(exists $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_id'});
    my $gnos_id = $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_id'};

    my $download = sprintf $GTDL_COMMAND, $options->{'gtdownload'},
                                          $options->{'gt_timeout'},
                                          $options->{'keys'}->{$repo},
                                          $repo,
                                          $gnos_id,
                                          $outbase;
    my $f_base = $outbase.'/'.$gnos_id;
    my $success = $f_base.'/.success.t';
    if(-e $success) {
      make_path($donor_base) unless(-e $donor_base);
      check_or_create_symlink($f_base, "$donor_base/$caller");
      next;
    }
    next if($options->{'symlinks'});

    my $out_file = "$f_base.out.log";
    my $err_file = "$f_base.err.log";
    $download .= " > $f_base.out.log";
    $download = "($download) >& $err_file";
    warn "Executing: $download\n";

    my $exit = system($download);

    return 0 if($exit && download_abandoned($err_file));

    die "A problem occured while executing: $download\n\tPlease check $f_base.err.log and proxy settings\n" if($exit);
    # as successful clean up logs
    unlink "$f_base.out.log";
    unlink "$f_base.err.log";

    # now build symlinks
    make_path($donor_base) unless(-e $donor_base);
    check_or_create_symlink($f_base, "$donor_base/$caller");

    # touch a success file in the output loc
    open my $S, '>', $success; close $S;
  }
  return 1;
}

sub download_abandoned {
  my $err_file = shift;
  my ($g_stdout, $g_stderr, $g_exit) = capture { system( sprintf q{grep -cF 'Inactivity timeout triggered after' %s}, $err_file); };
  chomp $g_stdout;
  if($g_exit == 0 && $g_stdout > 0) {
    warn "Abandoned download due to inactivity, skipping\n";
    return 1;
  }
  # a slightly different error
  ($g_stdout, $g_stderr, $g_exit) = capture { system( sprintf q{grep -cF 'Timeout was reached' %s}, $err_file); };
  chomp $g_stdout;
  if($g_exit == 0 && $g_stdout > 0) {
    warn "Abandoned download due gto timeout - could be problem with GNOS, skipping\n";
    return 1;
  }
  return 0;
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
  my %project_dist;
  my @filter_keys = keys %{$options->{'FILTERS'}};
  my $command = sprintf 'gunzip -c %s', $gz_manifest;
  my ($pid, $z);
  $pid = open $z, q{-|}, $command or croak 'Could not fork: '.$!;
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

    if(exists $options->{'COMPOSITE_FILTERS'}->{'multi_tumour'} && $donor->{'flags'}->{'all_tumor_specimen_aliquot_counts'} == 1) {
      warn sprintf "Donor: %s only has 1 tumour sample, requested 'COMPOSITE_FILTER.multi_tumour' is set\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }

    if(exists $options->{'COMPOSITE_FILTERS'}->{'manual_donor_blacklist'} && exists $options->{'COMPOSITE_FILTERS'}->{'manual_donor_blacklist'}->{ $donor->{'donor_unique_id'} }) {
      warn sprintf "Donor: %s has been defined in a filter 'COMPOSITE_FILTERS.manual_donor_blacklist'\n", $donor->{'donor_unique_id'} if($options->{'debug'});
      next;
    }


    if($options->{'analysis'} eq 'CALLS') {
      my @callers = @{$donor->{'flags'}->{'variant_calling_performed'}};
      my $keep = 0;
      for my $caller(@callers) {
        next if(exists  $options->{'COMPOSITE_FILTERS'}->{'caller'} && index(','.$options->{'COMPOSITE_FILTERS'}->{'caller'}.',', ','.$caller.',') == -1);
        next unless($donor->{'flags'}->{'is_'.$caller.'_variant_calling_performed'});
        next unless(exists $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'gnos_id'});
        my $size = data_size_calls_gb($options, $donor, $caller);
        $keep = 1;
        last;
      }
      unless($keep) {
        warn sprintf "Donor: %s has no variant calling available\n", $donor->{'donor_unique_id'} if($options->{'debug'});
        next;
      }

      $keep = 0;
      for my $caller(@callers) {
        next if($options->{'COMPOSITE_FILTERS'}->{'jamboree_approved'} == 1 && !$donor->{'flags'}->{'is_'.$caller.'_vcf_in_jamboree'});
        if(exists $options->{'COMPOSITE_FILTERS'}->{$caller.'_version'}) {
          my $result_version = $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'workflow_details'}->{'variant_workflow_version'};
          unless($result_version eq $options->{'COMPOSITE_FILTERS'}->{$caller.'_version'}) {
            warn sprintf "\tDonor: %s calls not from %s version %s\n", $donor->{'donor_unique_id'}, $caller, $options->{'COMPOSITE_FILTERS'}->{$caller.'_version'} if($options->{'debug'});
            next;
          }
          $keep++;
        }
        else {
          $keep++;
        }
        last;
      }
      unless($keep) {
        warn sprintf "Donor: %s has no variant calling available\n", $donor->{'donor_unique_id'} if($options->{'debug'});
        next;
      }

    }
    else {
      my $size = data_size_alignments_gb($options, $donor);
      if(exists $options->{'COMPOSITE_FILTERS'}->{'max_dataset_GB'} && $size > $options->{'COMPOSITE_FILTERS'}->{'max_dataset_GB'}) {
        warn sprintf "Donor: %s alignments exceed specified size limit of %s GB\n", $donor->{'donor_unique_id'}, $options->{'COMPOSITE_FILTERS'}->{'max_dataset_GB'} if($options->{'debug'});
        next;
      }
    }

    push @to_process, $donor;
    $project_dist{$donor->{'dcc_project_code'}}{'count'}++;
    $project_dist{$donor->{'dcc_project_code'}}{'total'} += $donor->{'_total_size'} || 0;
  }
  close $z;
  my $retained = scalar @to_process;
  print sprintf "Retained donors: %s\nRejected donors: %s\n", $retained, $total - $retained;
  print "Project Distribution\n";
  for(sort keys %project_dist) {
    my $avg = (int $project_dist{$_}{'total'} / $project_dist{$_}{'count'});
    $avg = $avg == 0 ? '?' : $avg+1;
    print "\t$_:\t$project_dist{$_}{count}\t(avg $avg GB)\n";
  }
  return \@to_process;
}

sub data_size_alignments_gb {
  my ($options, $donor) = @_;
  my $total_size = $donor->{'normal_alignment_status'}->{'aligned_bam'}->{'bam_file_size'};
  for my $tumour(@{$donor->{'tumor_alignment_status'}}) {
    $total_size += $tumour->{'aligned_bam'}->{'bam_file_size'};
  }
  $donor->{'_total_size'} = $total_size/1024/1024/1024;
  return $donor->{'_total_size'};
}

sub data_size_calls_gb {
  my ($options, $donor, $caller) = @_;
  my $total_size = 0;
  if(exists $donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'files'}) {
    for my $file(@{$donor->{'variant_calling_results'}->{$caller.'_variant_calling'}->{'files'}}) {
      $total_size += $file->{'file_size'};
    }
  }
  $donor->{'_total_size'} += $total_size/1024/1024/1024;
  return $donor->{'_total_size'};
}

sub manifest {
  my $options = shift;
  my $url = $options->{'url'} ? $options->{'url'} : $DEFAULT_URL;
  my $ff = File::Fetch->new(uri => $url, tempdir_root => File::Spec->tmpdir);
  my $listing;
  my $where = $ff->fetch( to => \$listing );
  my ($file) = $listing =~ m/(donor_p_[[:digit:]]+[.]jsonl[.]gz)/xms;
  my $to_get = "$url/$file";
  make_path($options->{'outdir'}) unless(-e $options->{'outdir'});
  $ff = File::Fetch->new(uri => $to_get);
  $where = $ff->fetch( to => $options->{'outdir'});
  my $dest = $options->{'outdir'}.'/'.$options->{'analysis'}.'_donor_p.jsonl.gz';
  move($where, $dest);
  return $dest;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'm|man' => \$opts{'m'},
		'i|info' => \$opts{'info'},
		's|symlinks' => \$opts{'symlinks'},
		'u|url=s' => \$opts{'url'},
		't|threads=i' => \$opts{'threads'},
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

  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
  PCAP::Cli::file_for_reading('config', $opts{'config'});
  pod2usage(-msg  => "\nERROR: Option '-a' is missing or invalid.\n", -verbose => 2,  -output => \*STDERR) unless(defined $opts{'analysis'} && first {$opts{'analysis'} eq $_} @ANALYSIS_TYPES);


	return \%opts;
}

__END__

=head1 NAME

gnos_pull.pl - retrieve/update analysis flow results on local systems.

=head1 SYNOPSIS

./gnos_pull.pl [-h] -u http://pancancer.info/gnos_metadata/latest/ -c gnos_pull.ini -o local_mirror/

  Required input:

    --analysis  (-a)  ALIGNMENTS or CALLS

    --outdir    (-o)  Where to save jsonl and resulting GNOS downloads

    --config    (-c)  Mapping of GNOS repos to permissions keys

  Other options:

    --symlinks  (-s)  Rebuild symlinks only.

    --threads   (-t)  Number of parallel GNOS retrievals.

    --url       (-u)  The base URL to retrieve jsonl file from
                        [http://pancancer.info/gnos_metadata/latest/]

    --info      (-i)  Just prints how many donor's will be included in pull and some stats.

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

=item i|info

Find out how much data you will potentially retrieve.  Can also be used to see how data is distributed.

This doesn't take into account data already downloaded.

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


