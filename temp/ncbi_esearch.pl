#!/usr/bin/perl -w
use strict;
use LWP::Simple;

my $htxt = <<HH;
perl $0 nuccore ON597627.1 [fasta] > ON597627.1.fasta

Sample commands:
  perl ncbi_esearch.pl nuccore ON597627.1 genbank > Crehmii_ctGenome.gb
  perl ncbi_esearch.pl nuccore ON597627.1 fasta   > Crehmii_ctGenome.fa
  perl ncbi_esearch.pl nuccore ON597627.1         > Crehmii_ctGenome.fa
  
HH

!@ARGV and die $htxt;

# Download PubMed records that are indexed in MeSH for both asthma and 
# leukotrienes and were also published in 2009.

my $db      = shift; # 'pubmed'
my $query   = shift; # 'asthma[mesh]+AND+leukotrienes[mesh]+AND+2009[pdat]';
my $rettype = shift; $rettype //= 'fasta';

#assemble the esearch URL
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

#post the esearch URL
my $output = get($url);

#parse WebEnv and QueryKey
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for ESearch-ESummary
#assemble the esummary URL
# $url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

#post the esummary URL
# my $docsums = get($url);
# print "$docsums";

### include this code for ESearch-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=$rettype&retmode=text";

#post the efetch URL
my $data = get($url);
print "$data";

