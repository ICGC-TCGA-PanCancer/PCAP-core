#!/use/bin/perl

use strict;
use LWP::Simple;
use XML::Simple qw(:strict);
use JSON;
use PCAP;
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);

my $options = &setup;
xml_to_bas($options);


sub xml_to_bas {
  my $options = shift;
  my $document = XMLin( get($options->{'uri'})
                      , ForceArray => 1
                      , KeyAttr => [],);

  my $bas_json = find_bas_json($document);
  json_to_bas_file($bas_json, $options->{'output'});
  return 1;
}

sub json_to_bas_file {
  my ($bas_json, $out_path) = @_;

  my $bas_data = decode_json $bas_json;
  my @metrics = @{$bas_data->{'qc_metrics'}};

  my @columns = bas_columns($metrics[0]);

  my $OUT;
  if(defined $out_path) {
    open $OUT, '>', $out_path;
  }
  else {
    $OUT = *STDOUT;
  }

  print $OUT join "\t", ('read_group_id', @columns);
  print $OUT "\n";
  for my $row(@metrics) {
    my @cols = ($row->{'read_group_id'});
    for my $col(@columns) {
      push @cols, $row->{'metrics'}->{$col};
    }
    print $OUT join "\t", @cols;
    print $OUT  "\n";
  }

  close $OUT;
}

sub bas_columns {
  my $first_record = shift;
  my @columns = sort keys $first_record->{'metrics'};
  return @columns;
}


sub find_bas_json {
  my $document = shift;
  for my $result(@{$document->{'Result'}}) {
    next unless(exists $result->{'analysis_xml'});
    for my $analysis_xml(@{$result->{'analysis_xml'}}) {
      next unless(exists $analysis_xml->{'ANALYSIS_SET'});
      for my $analysis_set(@{$analysis_xml->{'ANALYSIS_SET'}}) {
        next unless(exists $analysis_set->{'ANALYSIS'});
        for my $analysis(@{$analysis_set->{'ANALYSIS'}}) {
          next unless(exists $analysis->{'ANALYSIS_ATTRIBUTES'});
          for my $analysis_attributes(@{$analysis->{'ANALYSIS_ATTRIBUTES'}}) {
            next unless(exists $analysis_attributes->{'ANALYSIS_ATTRIBUTE'});
            for my $analysis_attribute(@{$analysis_attributes->{'ANALYSIS_ATTRIBUTE'}}) {
              next unless($analysis_attribute->{'TAG'}->[0] eq 'qc_metrics');
              return $analysis_attribute->{'VALUE'}->[0];
            }
          }
        }
      }
    }
  }
}

sub setup{
  my %opts;
  my @random_args;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'd|uri=s' => \$opts{'uri'},
              'o|output=s' => \$opts{'output'},
              '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  my $version = PCAP->VERSION;

  if(defined $opts{'v'}){
    print "Version: $version\n";
    exit;
  }

  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'m'});

  pod2usage(-message  => "\nERROR: unrecognised commandline arguments: ".join(', ',@random_args).".\n", -verbose => 1,  -output => \*STDERR) if(scalar @random_args) ;
  pod2usage(-message  => "\nERROR: d|uri must be defined.\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{'uri'});

  return \%opts;
}

__END__

=head1 NAME

xml_to_bas.pl - Generates a file containing read statistics for a given XML analysisFull URI.

=head1 SYNOPSIS

bam_stats.pl [options] [file...]

  Required parameters:
    -uri    -d    Same URI used by gtdownload
    -output -o    Name for output file. Defaults to STDOUT.

  Other:
    -help     -h   Brief help message.
    -man      -m   Full documentation.
    -version  -v   Prints the version number.

    xml_to_bas.pl -d https://gtrepo-ebi.annailabs.com/cghub/metadata/analysisFull/4e183691-ba1f-4103-a517-948f363928b8 -o file.bam.bas

=head1 OPTIONS

=over 8

=item B<-uri>

Which BAS data to download and convert.

=item B<-output>

File path to output data. If this option is omitted the script will attempt to write to
STDOUT.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-version>

Prints the version number and exits.

=back

=cut
