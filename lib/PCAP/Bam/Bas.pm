package PCAP::Bam::Bas;

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014 ICGC PanCancer Project
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
use English qw( -no_match_vars );
use warnings FATAL=>'all';
use autodie qw( :all );
use Carp qw(croak carp);

sub new {
  my ($class, $bas) = @_;
  my $self = { };
  bless $self, $class;
  $self->_init($bas);
  return $self;
}

sub _init {
  my ($self, $bas) = @_;
  croak "No bas file defined" if(!defined $bas);
  die "*.bas file: $bas does not exist" unless(-e $bas);
  die "*.bas file: $bas is empty" unless(-s $bas);
  open my $IN, '<', $bas;
  $self->bas_keys($IN);
  $self->_import_data($IN);
  close $IN;
  return 1;
}

sub _import_data {
  my ($self, $fh) = @_;
  while(my $line = <$fh>) {
    chomp $line;
    my @bits = split /\t/, $line;
    my %rg;
    for my $key(@{$self->bas_keys}) {
      $rg{$key} = $bits[$self->{'key_pos_map'}->{$key}];
    }
    $self->{'_data'}->{$rg{'readgroup'}} = \%rg;
  }
  return 1;
}

sub bas_keys {
  my ($self, $key_fh) = @_;
  croak "bas_keys should only be initialised once\n" if(exists $self->{'keys'} && defined $key_fh);
  if(defined $key_fh) {
    my $line = <$key_fh>;
    chomp $line;
    my @head = split /\t/, $line;
    my %key_pos_map;
    my $pos=0;
    for my $key(@head) {
      $key_pos_map{$key} = $pos++;
    }
    $self->{'keys'} = \@head;
    $self->{'key_pos_map'} = \%key_pos_map;
  }
  return $self->{'keys'};
}

sub read_groups {
  return (sort keys shift->{'_data'});
}

sub get {
  my ($self, $rg, $key) = @_;
  die qq{Readgroup '$rg' does not exist\n} unless(exists $self->{'_data'}->{$rg});
  return exists $self->{'_data'}->{$rg}->{$key} ? $self->{'_data'}->{$rg}->{$key} : undef;
}

1;

__END__

=head1 PCAP::Bam::Bas

Convenience class for accessing data in a *.bas file.

=head2 METHODS

=over 2

=item new

Construct an access object for BAM statistics file.

 my $bas_ob = PCAP::Bam::Bas->new($bas);

=item bas_keys

Returns the list of available keys for this BAS file.

=item read_groups

Returns sorted list of read-groups found in this BAS file.

=item get

Retrieve a value by its readgroup and key:

 $bas->($rg, 'median_insert_size');

NOTE: Returns undef if a key is not available.

=back
