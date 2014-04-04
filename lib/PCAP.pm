package PCAP;

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


use strict;
use Const::Fast qw(const);

our $VERSION = '0.2.99';

const my $LICENSE =>
"#################
# PCAP version %s, Copyright (C) 2014 ICGC/TCGA Pan-Cancer Analysis Project
# PCAP comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

const my $DEFAULT_PATH => 'biobambam,samtools,bwa';
const my %UPGRADE_PATH => ( '0.1.0' => 'biobambam,samtools,bwa',
                            '0.1.1' => 'biobambam,bwa',
                            '0.1.2' => 'biobambam',
                            '0.2.0' => 'biobambam',
                            '0.3.0' => '',
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
