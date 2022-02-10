#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 ipr.SearchResults-transposases.json > ipr.SearchResults-transposases.json.tbl\n";

print STDOUT join("\t", qw/id source source_database name description/)."\n";
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @lines = &parse_iprJson($_);
  for my $a1 (@lines) {
    for my $k1 (keys %$a1) {
      $a1->{$k1} //= "NA";
    }
    print STDOUT join("\t", @{$a1}{qw/id source source_database name description/})."\n";
  }
}

sub parse_iprJson {
  my ($txt1) = @_;
  $txt1 =~ s!^\s*\[(.*)\]\s*$!$1! or die "Err 1:$txt1\n";
  my @back;
  while ($txt1 =~ s!^\s*\{ "id":"(\S+?)", "source":"([^"]+)", "fields":\{ "description":\[(?:"(.*?)")?\], "name":\["([^"]+)"\], "source_database":\["([^"]+)"\] \} \s* \}\s*,*!!x) {
    my %h;
    @h{qw/id source description name source_database/} = ($1, $2, $3, $4, $5);
    push(@back, \%h);
  }
  $txt1 =~ m!^\s*$! or die "Err 2: |$txt1|\n";
  return(@back);
}# parse_iprJson()

# {
#"id":"PTHR22955",
#"source":"interpro7_family",
#"fields":{
#  "description":[],
#  "name":["RETROTRANSPOSON"],
#  "source_database":["PANTHER"]
#}
#}

