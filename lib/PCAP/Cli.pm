package PCAP::Cli;

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the IGCG/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014 IGCG PanCancer Project
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


use PCAP;
our $VERSION = PCAP->VERSION;

use strict;
use autodie qw(:all);
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use File::Path qw (make_path);
use List::Util qw(first);

sub file_for_reading {
  my ($opt_name, $opt_val) = @_;
  die "Option '$opt_name' has not been defined.\n" unless(defined $opt_val);
  die "Option '$opt_name' ($opt_val) should be an existing file.\n" unless(-e $opt_val);
  die "Option '$opt_name' ($opt_val) is a directory.\n" if(-d $opt_val);
  die "Option '$opt_name' ($opt_val) should be file with non-zero size.\n" unless(-s $opt_val);
  return $opt_val;
}

sub out_dir_check {
  my ($opt_name, $opt_val) = @_;
  die "Option '$opt_name' has not been defined.\n" unless(defined $opt_val);

  if(-e $opt_val) {
    if(-d $opt_val) {
      # check writable
      die "Option '$opt_name' points to an existing WRITE PROTECTED directory: $opt_val\n." unless(-w $opt_val);
    }
    else {
      die "Option '$opt_name' points to an existing entity (not a directory): $opt_val\n.";
    }
  }
  else {
    make_path($opt_val, {error => \my $err});
    if(@{$err}) {
      for my $diag (@$err) {
        my ($file, $message) = %$diag;
        # uncoverable branch false
        $message .= sprintf ' (%s)', $file if($file ne {});
        die "$message.\n";
      }
    }
  }
  return $opt_val;
}

sub valid_process {
  my ($opt_name, $opt_val, $valid_procs) = @_;
  die "Option '$opt_name' is not an expected process type: $opt_val\n." unless(first {$opt_val eq $_} @{$valid_procs});
  return $opt_val;
}

sub valid_index_by_factor {
  my ($opt_name, $opt_val, $base, $proc_factor) = @_;
  $proc_factor ||= 1;
  my $max = $base * $proc_factor;
  die "Option '$opt_name' needs to be between 1 and $max: $opt_val\n." unless($opt_val > 0 && $opt_val <= $max);
  return $opt_val;
}

sub opt_requires_opts {
  my ($opt_name, $options, $req_ops) = @_;
  my $count = 0;
  for(@{$req_ops}) {
    $count++ if(exists $options->{$_});
  }
  die "Option '$opt_name' requires these options: ".(join q{,},@{$req_ops})."\n" if($count != scalar @{$req_ops});
  return 1;
}

1;

__END__

=head1 NAME

PCAP::Cli - Collection of functions for checking inputs.

=head2 Methods

=over 4

=item file_for_reading

Checks input is an existing file with size.

  my $checked = PCAP::Cli::file_for_reading($option_name, $option_value);

The name of the option is used to aid in generating informative messages to the user.

=item out_dir_check

Checks if provided path is present and writable, if not will attempt to create.

  my $checked = PCAP::Cli::out_dir_check($option_name, $option_value);

The name of the option is used to aid in generating informative messages to the user.

=back
