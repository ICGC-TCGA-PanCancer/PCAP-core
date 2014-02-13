package PCAP::Bam;

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
use autodie qw(:all);
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Const::Fast qw(const);
use File::Spec;
use File::Which qw(which);
use Bio::DB::Sam;
use Carp qw(croak);
use List::Util qw(first);

use PCAP::Threaded;

const my $BAMCOLLATE => q{(%s colsbs=268435456 collate=1 reset=1 exclude=SECONDARY,QCFAIL,SUPPLEMENTARY T=%s filename=%s level=1 > %s)};
const my $BAMBAM_DUP => q{ O=%s M=%s tmpfile=%s markthreads=%s rewritebam=1 rewritebamlevel=1 index=1 md5=1};

sub new {
  my ($class, $bam) = @_;
  my $self = {};
  if(defined $bam) {
    $self->{'bam'} = $bam;
    $self->{'md5'} = $bam.'.md5' if(-e $bam.'.md5');
  }
  bless $self, $class;
  return $self;
}

sub rg_line_for_output {
  my $bam = shift;
  my $sam = sam_ob($bam);
  my $header = $sam->header->text;
  my $rg_line;
  while($header =~ m/^(\@RG\t[^\n]+)/xmsg) {
    my $new_rg = $1;
    die "BAM file appears to contain data for multiple readgroups, not supported: \n\n$header\n" if(defined $rg_line);
    $rg_line = $new_rg;
    $rg_line =~ s/\t/\\t/g;
  }
  return ($rg_line, $sam); # also return the SAM object
}

sub bam_to_grouped_bam {
  # uncoverable subroutine
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  my $tmp = $options->{'tmp'};
  # uncoverable branch true
  # uncoverable branch false
  return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
  my $inputs = $options->{'meta_set'};
  my $bamcollate =  which('bamcollate2') || die "Unable to find 'bamcollate2' in path";
  my $command = sprintf $BAMCOLLATE, $bamcollate, File::Spec->catfile($tmp, "collate.$index"), $inputs->[$index-1]->in, $inputs->[$index-1]->tstub.'.bam';
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub merge_and_mark_dup {
  # uncoverable subroutine
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my $marked = File::Spec->catdir($options->{'outdir'}, $options->{'sample'});
  my $met = "$marked.met";
  $marked .= '.bam';
  my $helper_threads = $options->{'threads'}-1;
  # uncoverable branch true
  # uncoverable branch false
  $helper_threads = 1 if($helper_threads < 1);
  my $command = which('bammarkduplicates') || die "Unable to find 'bammarkduplicates' in path";
  $command .= sprintf $BAMBAM_DUP,  $marked,
                                    $met,
                                    File::Spec->catfile($tmp, 'biormdup'),
                                    $helper_threads;
  for(@{$options->{'meta_set'}}) {
    $command .= ' I='.$_->tstub.'_sorted.bam';
  }
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
  return $marked;
}

sub sample_name {
  my $bam = shift;
  my $sam = sam_ob($bam);
  my $header = $sam->header->text;
  my $sample;
  while($header =~ m/\tSM:([^\t\n]+)/xmsg) {
    my $new_sample = $1;
    die "BAM file appears to contain data for multiple samples, not supported: \n\n$header\n" if(defined $sample && $sample ne $new_sample);
    $sample = $new_sample;
  }
  return ($sample, $sam); # also return the SAM object
}

sub read_group_info {
  my ($self, $required_tags) = @_;

  # sort for convention rather than any purpose
  # de-reference for ease
  my @expected_tags;
  @expected_tags = sort @{$required_tags} if(defined $required_tags);

  if(exists $self->{'rgs'}) {
    $self->check_for_tags(\@expected_tags) if(defined $required_tags);
    return $self->{'rgs'};
  }

  my @rgs;
  for my $line (@{$self->sam_header}) {
    next unless($line =~ m/^\@RG/);
    my @elements = split /\t/, $line;
    shift @elements; # drop @RG entry
    my %rg;
    for my $element(@elements) {
      # Don't think it's possible to generate header with blank fields (\t+)
      # but lets be parnoid
      next if($element eq q{});

      my ($tag, $value) = $element =~ m/^([^\:]+)\:(.+)$/;
      die "ERROR: Malformed RG tag/value: $element\n" unless(defined $tag && defined $value);
      $rg{$tag} = $value;
    }
    push @rgs, \%rg;
  }
  die "ERROR: This BAM has no readgroups: $self->{bam}" if(scalar @rgs == 0);
  $self->{'rgs'} = \@rgs;
  $self->check_for_tags(\@expected_tags) if(defined $required_tags);
  return $self->{'rgs'};
}

sub comments {
  my $self = shift;
  my @comments;
  for my $line(@{$self->sam_header}) {
    next unless($line =~ m/^\@CO\t(.*)/);
    push @comments, $1;
  }
  return \@comments;
}

sub sam_header {
  my $self = shift;
  return $self->{'header'} if(exists $self->{'header'});
  my $sam = sam_ob($self->{'bam'});
  my @header_lines = split /\n/, $sam->header->text;
  $self->{'header'} = \@header_lines;
  return $self->{'header'};
}

sub single_rg_value {
  my ($self, $tag) = @_;
  my $rgs = $self->read_group_info([$tag]);
  croak "ERROR: This BAM includes multiple readgroups: $self->{bam}" if(scalar @{$rgs} > 1);
  return $rgs->[0]->{$tag};
}

sub check_for_tags {
  my ($self, $sorted_tags) = @_;
  my @errors;
  for my $rg(@{$self->{'rgs'}}) {
    my @tags_found = keys %{$rg};
    for my $i(0..(scalar @{$sorted_tags})-1) {
      unless(first { $sorted_tags->[$i] eq $_ } @tags_found) {
        my @pairs;
        push @pairs, "$_:$rg->{$_}" for(@tags_found);
        push @errors, "ERROR: $sorted_tags->[$i] not found in RG of $self->{bam}\n\t".(join q{,}, @pairs);
      }
    }
  }
  die (join "\n", @errors) if(scalar @errors > 0);
  return 1;
}

sub sam_ob {
  my $bam = shift;
  my $sam;
  if(ref $bam eq 'Bio::DB::Sam') {
    $sam = $bam;
  }
  elsif(-e $bam) {
    $sam = Bio::DB::Sam->new(-bam => $bam);
  }
  else {
    my $caller = (caller(1))[3];
    croak "$caller requires either a BAM file or Bio::DB::Sam object.\n";
  }
  return $sam;
}

sub check_paired {
  my $self = shift;
  my $sam = sam_ob($self->{'bam'});
  my $bam = $sam->bam;
  $bam->header;
  my $read = $bam->read1;
  die "ERROR: Input BAMs should be for paired end sequencing: $self->{bam}\n" unless(1 & $read->flag);
  return 1;
}

1;

__END__

=head1 NAME

PCAP::Bam - Methods that process BAM files

=head2 Object Methods

=over 4

=item new

  my $bam = PCAP::Bam->new($bam_file);

Generate a BAM object to allow access to object based functions.  Primarily added to allow efficient
access to header information without re-parsing once populated.

=item read_group_info

  my $rg_data = $bam->read_group_info($required_tags);

Returns an arrayref (one entry per readgroup) of hashes where the keys are the RG tags found in the @RG header line.
Providing an arrayref of required_tags causes an error to be raised if any are not found:

  my $rg_data = $bam->read_group_info([qw(CN LB SN)]);

Resulting data structure:

  [ { CN => 'SI',
      LB => 'LIBRARY_ID',
      SM => 'Sample_name},
    { CN => ...
    },
    ...];

=item single_rg_value

  my $sample_name = $bam->single_rg_value('SM');

Gets the value of a single RG tag.  Errors if multiple readgroups detected.

=item check_for_tags

  $bam->check_for_tags([list of required tags]);

Check that all readgroups in this BAM have the specified tags.
Errors if any are not detected.

=item check_paired

  $bam->check_paired;

Will error if BAM file doesn't contain paired reads.

=item comments

  my @comments = @{$bam->comments};

Returns an array ref of the value of each comment line.
'@CO\t' is pre-stripped.

=item sam_header

  my @header_lines = @{$bam->sam_header};

Returns array ref of all header lines.  One entry perl line.

=back

=head2 Non Object Methods

=over 4

=item bam_to_grouped_bam

  PCAP::Bam::bam_to_grouped_bam($index, $options);

Convert a BAM for a single readgroup to a readname grouped BAM.

  index   - Which element of the 'bams' array reference found in options should be processed
              NOTE: index is origin 1

  options - Hashref, requires the following entries:

          -tmp  : working/output directory depending on application
          -bams : array reference to a B<sorted> list of BAM files (to allow resume function)

On successful completion of a run a file is created in tmp/progress/bamcollate.$index.  If the
program is re-run without clearing of this progress folder the execution of the command is skipped.

=item merge_and_mark_dup

  PCAP::Bam::merge_and_mark_dup($options);

Takes a list of sorted BAM files, marks duplicates and produces single output BAM in one pass.

  options - Hashref, requires the following entries:

          -tmp      : working/output directory depending on application
          -bams     : array reference to a B<sorted> list of BAM files (for consistent results)
          -outdir   : output files written to this folder
          -sample   : prefix of output file
          -threads  : Total threads available to process

  returns - path_to_marked_bam

Resulting data is of the form:

  $outdir/$sample.bam
  $outdir/$sample.bam.bai
  $outdir/$sample.bam.md5
  $outdir/$sample.met

=item sample_name

Takes BAM or Bio::DB::Sam object as input and returns the sample name found in the header.

The SAM object is also returned should it be useful for other calls

=item rg_line_for_output

Takes BAM or Bio::DB::Sam object as input and returns the string representation for the RG line.
Intended for use when adding RG to BWA MEM output and is only useful in single RG BAMs

The SAM object is also returned should it be useful for other calls

=item sam_ob

  my $sam_ob = sam_ob('file.bam');

Generate a Bio::DB::Sam object from the provided BAM file.

=back
