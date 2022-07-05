#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and -t and die "perl $0 in.bn6 > in.bn6.1query1line\n";

my @bn6_colName = qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle/;
my @col_want = (1,2,3,10,11,17,16,18);
my $topN = 5;
my @txt_want = (@bn6_colName[@col_want]) x $topN;
print STDOUT join("\t", "qseqid", @txt_want)."\n";



my (%h, @qIDs);
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  defined $h{$ta[0]} or push(@qIDs, $ta[0]);
  push(@{$h{$ta[0]}}, [$ta[11], [@ta[@col_want]]]);
}
for my $qid (@qIDs) {
  my @out_line = ($qid);
  @{$h{$qid}} = sort { $b->[0] <=> $a->[0] } @{$h{$qid}};
  my $i=0;
  for my $t1 (@{$h{$qid}}) {
    $i < $topN or last;
    push(@out_line, @{$t1->[1]});
    $i++;
  }
  print STDOUT join("\t", @out_line)."\n";
}

