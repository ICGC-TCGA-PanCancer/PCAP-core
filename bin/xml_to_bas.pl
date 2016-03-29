#!/use/bin/perl

use strict;
use Capture::Tiny qw(capture);
use XML::Simple qw(:strict);
use JSON;
use PCAP;
use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Bio::DB::HTS;
use FindBin qw($Bin);
use lib "$Bin/../lib";

my $options = &setup;
xml_to_bas($options);

sub get {
  my $uri = shift;
   if ($uri =~ /http:\/\// || $uri =~ /https:\/\//) {
     my ($raw_xml, $e_str, $e_code) = capture { system('curl -ksS '.$uri); };
     die "ERROR: $e_str\n" if($e_code);
     return $raw_xml;
   } else {
     my ($raw_xml, $e_str, $e_code) = capture { system('cat '.$uri); };
     die "ERROR: $e_str\n" if($e_code);
     return $raw_xml;
   }
}

sub xml_to_bas {
  my $options = shift;
  my $input_path = $options->{'uri'};
   if (defined $options->{'local-path'}) {
     $input_path = $options->{'local-path'};
   }
  my $document = XMLin( get($input_path)
                      , ForceArray => 1
                      , KeyAttr => [],);

  my $bas_json = find_bas_json($document);
  json_to_bas_file($bas_json, $options);
  return 1;
}

sub json_to_bas_file {
  my ($bas_json, $options) = @_;

  my $bas_data = decode_json $bas_json;
  validate_bas($bas_data, $options);
  my @metrics = @{$bas_data->{'qc_metrics'}};

  my @columns = bas_columns($metrics[0]);

  my $OUT;
  if(defined $options->{'output'}) {
    open $OUT, '>', $options->{'output'};
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
  my @columns = sort keys %{$first_record->{'metrics'}};
  return @columns;
}

sub validate_bas {
  my ($bas_data, $options) = @_;
  # first look for read_group_id clash
  my %ids;
  my $clash = clash_check($bas_data);
  if($clash) {
    unless($options->{'bam'}) {
      die "ERROR: multiple metric entries with the same read_group_id. May be recoverable with '-b' defined.\n";
    }
    else {
      warn "WARNING: multiple metric entries with the same read_group_id attempting to compensate...\n";
      correct_clash($bas_data, $options->{'bam'});
    }
  }
}

sub correct_clash {
  my ($bas_data, $bamfile) = @_;
  my $bam = Bio::DB::HTS->new(-bam => $options->{'bam'});
  my %rg_by_pu;
  foreach my $hl(split /\n/, $bam->header->text) {
    next unless($hl =~ m/^\@RG/);
    my ($pu) = $hl =~ m/\tPU:([^\t]+)/;
    my ($id) = $hl =~ m/\tID:([^\t]+)/;
    die "ERROR: Unable to recover read_group_id clash using PU field.\n" if(exists $rg_by_pu{$pu});
    $rg_by_pu{$pu} = $id;
  }

  my @metrics = @{$bas_data->{'qc_metrics'}};
  for my $row(@metrics) {
    $row->{'read_group_id'} = $rg_by_pu{ $row->{'metrics'}->{'platform_unit'} };
    $row->{'metrics'}->{'readgroup'} = $rg_by_pu{ $row->{'metrics'}->{'platform_unit'} };
  }
  die "ERROR: Unable to recover read_group_id clash using PU field.\n" if(clash_check($bas_data));
  return 1;
}

sub clash_check {
  my $bas_data = shift;
  my @metrics = @{$bas_data->{'qc_metrics'}};
  my %ids;
  my $clash = 0;
  for my $row(@metrics) {
    if(exists $ids{$row->{'read_group_id'}}) {
      $clash++;
    }
    else {
      $ids{$row->{'read_group_id'}} = 1;
    }
  }
  return $clash;
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
              'l|local-path=s' => \$opts{'local-path'},
              'b|bam=s' => \$opts{'bam'},
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

  pod2usage(-message  => "\nERROR: unrecognised commandline arguments: ".join(', ',@random_args).".\n", -verbose => 1,  -output => \*STDERR) if(scalar @random_args);
  pod2usage(-message  => "\nERROR: either d|uri or l|local-path must be defined.\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{'uri'} || defined $opts{'local-path'});

  return \%opts;
}

__END__

=head1 NAME

xml_to_bas.pl - Generates a file containing read statistics for a given XML analysisFull URI.

=head1 SYNOPSIS

xml_to_bas.pl [options]

  Required parameters:
    -uri        -d    Same URI used by gtdownload
    -output     -o    Name for output file. Defaults to STDOUT.

  Optional parameters:
    -bam       -b  BAM|CRAM file this data relates to
                   - checks retrieved data correlates with expected BAM
                   - additionally can correct read_group_id if other fields correlate when
                     clashes occur.
    -local-path -l Local file path for XML from GNOS.

  Other:
    -help     -h   Brief help message.
    -man      -m   Full documentation.
    -version  -v   Prints the version number.

  Example:
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
