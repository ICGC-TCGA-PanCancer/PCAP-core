#!/usr/bin/env perl

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
use warnings FATAL => 'all';
use autodie qw(:all);
use English qw( -no_match_vars );

use Proc::ProcessTable;
use List::Util qw(sum0 first);
use Const::Fast qw(const);

const my @OK_STATES => qw(run sleep);

die "ERROR: Requires a command to monitor\n" if(scalar @ARGV == 0);

my $started = time;
my $forked_id = fork();
if($forked_id == 0) {
  system(join q{ }, @ARGV);
}
else {
  my %processes;
  my ($rss_max, $virt_max) = (0,0);

  my $forked_running = 1;
  my $pt = new Proc::ProcessTable();
  while($forked_running == 1) {
    $forked_running = 0; # need to set if still running
    my (@rss_mem, @virt_mem);
    for my $p(@{$pt->table}) {
      # I only care about this process group
      next unless ($p->pgrp == $PROCESS_ID);

      # keep running until child ends
      $forked_running = 1 if($p->pid == $forked_id && $p->state ne 'defunct');

      push @rss_mem, $p->rss;
      push @virt_mem, $p->size;

      # capture time variables for each all processes in group
      $processes{$p->pid} = { 'stime' => $p->stime,
                              'utime' => $p->utime,};
    }
    # calculate total memory NOW, record if an increase
    my $rss_new = sum0 @rss_mem;
    $rss_max = $rss_new if($rss_new > $rss_max);
    my $virt_new = sum0 @virt_mem;
    $virt_max = $virt_new if($virt_new > $virt_max);

    sleep 1;
  }

  my ($stime, $utime);
  for my $p(keys %processes) {
    $stime += $processes{$p}{'stime'};
    $utime += $processes{$p}{'utime'};
  }

  my $output_format = <<OUTFORMAT;

####### JOB OUTPUT ABOVE THIS LINE #######

JOB SUMMARY
---------------------
Wall sec.    : %d
Time sec.    : %d
User sec.    : %d
System sec.  : %d
Peak Res GB  : %.2f
Peak Virt GB : %.2f

Time, User and System are a total of the time taken by all processes in each state.
Memory is peak based on sampling of all memory in use at 1 second intervals.
OUTFORMAT

  print sprintf $output_format, time - $started,
                                ($utime+$stime) /1000000,
                                $utime /1000000,
                                $stime /1000000,
                                $rss_max /1024/1024/1024,
                                $virt_max /1024/1024/1024;
}
