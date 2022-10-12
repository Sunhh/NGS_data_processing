#!/usr/bin/perl
use strict;
use warnings;

scalar(@ARGV) == 2 or die "perl $0 input/grp_tag comb.grp1 > comb.grp1.fmt\n";

my $f1_spec2tag = shift;
my (%tag2spec, @spec_a);
{
  my %spec_h;
  open F1,'<',"$f1_spec2tag" or die;
  while (<F1>) {
    chomp;
    my @ta=split(/\t/, $_);
    $tag2spec{$ta[1]} = $ta[0];
    $spec_h{$ta[0]} = (defined $ta[2]) ? $ta[2] : $.;
  }
  close F1;
  @spec_a = sort { $spec_h{$a} <=> $spec_h{$b} } keys %spec_h;
}

print STDOUT join("\t", "Grp_name", "size", (@spec_a), (map { "${_}_geneIDs" } @spec_a))."\n";
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my (%l_cont);
  for my $tb (@ta[2..$#ta]) {
    # $tb =~ s!^([^\s:]{1,5}):!! or die "|$tb|\n";
    $tb =~ s!^([^\s:]+):!! or die "Bad name |$tb|\n";
    my $sp_tag = $1;
    defined $tag2spec{$sp_tag} or die "bad tag [$sp_tag] in $tb\n";
    my $sp_txt = $tag2spec{$sp_tag};
    push(@{$l_cont{$sp_txt}}, $tb);
  }
  my (@spSize, @spCont);
  for my $sp (@spec_a) {
    $l_cont{$sp} //= [];
    push(@spSize, scalar(@{$l_cont{$sp}}));
    push(@spCont, join(";", @{$l_cont{$sp}}));
  }
  print STDOUT join("\t", $ta[0], $ta[1], @spSize, @spCont)."\n";
}

