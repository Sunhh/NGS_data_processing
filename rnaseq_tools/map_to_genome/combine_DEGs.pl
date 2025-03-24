#!/usr/bin/perl
use strict;
use warnings;

# Check input
die "Usage: $0 meta_file > combined_output.tsv\n" unless @ARGV == 1;
my $meta_file = $ARGV[0];

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

# Read each DEG file
foreach my $file (@files) {
    open IN, "<", $file or die "Cannot open $file: $!";
    my $header = <IN>; # Skip header line
    while (<IN>) {
        chomp;
        my ($gene, @vals) = split /\t/;
        push @gene_ids, $gene if !$header_found;
        push @{ $data{$gene} }, @vals;
    }
    $header_found = 1;
    close IN;
}

# Print headers
print "GeneID";
foreach my $label (@labels) {
    print "\t$label\t$label\t$label\t$label";
}
print "\n";

# Print sub-header
print "GeneID";
foreach my $f (@files) {
    print "\tMeanTPM_Baseline\tMeanTPM_Treatment\tFDR\tLog2FoldChange";
}
print "\n";

# Print data
foreach my $gene (@gene_ids) {
    print $gene;
    print "\t", join("\t", @{ $data{$gene} });
    print "\n";
}
