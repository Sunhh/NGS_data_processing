#!/usr/bin/perl
use strict;
use warnings;

# Check input
die "Usage: $0 meta_file  [columns:2-3,5] > combined_output.tsv\n" unless @ARGV >= 1;

my $meta_file = $ARGV[0];
my @coln;
scalar(@ARGV) > 1 and @coln = &parseCol($ARGV[1]);
for (@coln) { $_--; $_ >= 0 or die "Cannot add 0-col\n"; }

# Store file paths and labels
my @files;
my @labels;

open META, "<", $meta_file or die "Cannot open $meta_file: $!";
while (<META>) {
    chomp;
    next if /^\s*(#|\s*$)/;
    my ($file, $label) = split /\s+/, $_;
    push @files, $file;
    push @labels, $label;
}
close META;

# Hash to store all data per gene
my %data;
my @gene_ids;
my $header_found = 0;
my ($first_name, @header_names);

# Read each DEG file
foreach my $file (@files) {
  open IN, "<", $file or die "Cannot open $file: $!";
  my $header = <IN>;
  chomp($header);
  ($first_name, @header_names) = split(/\t/, $header);
  if (scalar(@coln) == 0) {
    @coln = (0 .. $#header_names);
  }
  @header_names = @header_names[@coln];
  while (<IN>) {
    chomp;
    my ($gene, @vals) = split /\t/;
    @vals = @vals[@coln];
    push @gene_ids, $gene if !$header_found;
    push @{ $data{$gene} }, @vals;
  }
  $header_found = 1;
  close IN;
}

# Print headers
print "$first_name";
foreach my $label (@labels) {
  print "\t$label" x scalar(@coln);
}
print "\n";

# Print sub-header
print "$first_name";
foreach my $f (@files) {
  print "\t".join("\t", @header_names);
}
print "\n";

# Print data
foreach my $gene (@gene_ids) {
    print $gene;
    print "\t", join("\t", @{ $data{$gene} });
    print "\n";
}

sub parseCol {
  my @cols = split(/,/, $_[0]);
  my @ncols;
  for my $tc (@cols) {
    $tc =~ s/(^\s+|\s+$)//g;
    if ($tc =~ m/^\d+$/) {
      push(@ncols, $tc);
    } elsif ($tc =~ m/^(\d+)\-(\d+)$/) {
      my ($s, $e) = ($1, $2);
      if ($s <= $e) {
        push(@ncols, ($s .. $e));
      }else{  
        push(@ncols, reverse($e .. $s));
      }
    } else {
      die("[Err]Unparsable column tag for [$tc]\n");
    }
  }
  return (@ncols);
}

