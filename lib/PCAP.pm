package PCAP;

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014-2016 ICGC PanCancer Project
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
use Const::Fast qw(const);
use base 'Exporter';
use FindBin qw($Bin);
use File::Which qw(which);

our $VERSION = '3.0.1';
our @EXPORT = qw($VERSION _which);

const my $LICENSE =>
"#################
# PCAP version %s, Copyright (C) 2014-2016 ICGC/TCGA Pan-Cancer Analysis Project
# PCAP comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

const my $DEFAULT_PATH => 'biobambam,samtools,bwa';
const my %UPGRADE_PATH => ( # all earlier versions need full upgrade
                            '1.2.0'  => 'biobambam,samtools,bwa', # if later versions have new versions then all preceding need that tool listing
                            '1.2.1'  => 'biobambam,samtools,bwa',
                            '1.2.2'  => 'biobambam,samtools,bwa',
                            '1.3.0'  => 'biobambam,samtools,bwa',
                            '1.4.0'  => 'biobambam,samtools,bwa',
                            '1.5.0'  => 'biobambam,samtools,bwa',
                            '1.5.1'  => 'biobambam,samtools,bwa',
                            '1.5.2'  => 'biobambam,samtools,bwa',
                            '1.5.3'  => 'biobambam,samtools,bwa',
                            '1.5.4'  => 'biobambam,samtools,bwa',
                            '1.6.0'  => 'biobambam,samtools,bwa',
                            '1.6.1'  => 'biobambam,samtools,bwa',
                            '1.6.2'  => 'biobambam,samtools,bwa',
                            '1.6.3'  => 'biobambam,samtools,bwa',
                            '1.7.0'  => 'biobambam,samtools,bwa',
                            '1.7.1'  => 'biobambam,samtools,bwa',
                            '1.8.0'  => 'biobambam,samtools',
                            '1.8.1'  => 'biobambam,samtools',
                            '1.8.2'  => 'biobambam,samtools',
                            '1.9.0'  => 'biobambam,samtools',
                            '1.9.1'  => 'biobambam,samtools',
                            '1.9.4'  => 'biobambam,samtools',
                            '1.10.0'  => 'biobambam,samtools',
                            '1.11.0'  => 'biobambam,samtools',
                            '1.11.1'  => 'biobambam,samtools',
                            '1.12.0'  => 'biobambam,samtools',
                            '1.12.1'  => 'biobambam,samtools',
                            '1.12.2'  => 'biobambam,samtools',
                            '1.12.3'  => 'biobambam,samtools',
                            '1.13.0'  => 'biobambam,samtools',
                            '1.13.1'  => 'biobambam,samtools',
                            '1.13.2'  => 'biobambam,samtools',
                            '1.14.0'  => 'biobambam,samtools',
                            '2.0.0'  => 'biobambam',
                            '2.0.1' => 'biobambam,samtools',
                            '2.1.0' => 'biobambam,samtools',
                            '2.1.2' => 'biobambam,samtools',
                            '2.1.3' => 'biobambam,samtools',
                            '2.2.0' => 'biobambam,samtools',
                            '2.2.1' => 'biobambam',
                            '2.3.0' => 'biobambam',
                            '2.4.0' => 'biobambam',
                            '2.5.0' => 'biobambam',
                            '3.0.0' => '',
                          );

sub license {
  return sprintf $LICENSE, $VERSION;
}

sub upgrade_path {
  my $installed_version = shift;
  return $DEFAULT_PATH if(!defined $installed_version);
  chomp $installed_version;
  return $DEFAULT_PATH if(!exists $UPGRADE_PATH{$installed_version});
  return $UPGRADE_PATH{$installed_version};
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
  return $path;
}

1;

__END__

=head1 NAME

PCAP - Base class to house version and generic functions.

=head2 Methods

=over 4

=item license

  my $brief_license = PCAP::licence;

Output the brief license text for use in help messages.

=item upgrade_path

  my $install_these = PCAP::upgrade_path('<current_version>');

Return the list of tools that should be installed by setup.sh when upgrading from a previous version.

=back
