#!/usr/bin/perl
use strict;
use warnings;
use SeqAlnSunhh;
use Getopt::Long;

my %opts;
GetOptions(\%opts,
  "nameSorted!",
);

-t and !@ARGV and die "samtools view in.bam | perl $0 [-nameSorted] > list.rdID_R12_rdLen_rdMismatch\n";

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) };
my %flag_P  = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,6=1;0=0')};

print STDOUT join("\t", qw/readName R12 readLength mismatch/)."\n";
if ($opts{'nameSorted'}) {
  my ($prevID, $prevLen, $prevMis) = ('', '', '');
  while (<>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[5] eq '*' and next;
    defined $flag_UN{$ta[1]} and next; # Read is unmapped.
    my $rdID = (defined $flag_P{$ta[1]}) ? "$ta[0]|R1" : "$ta[0]|R2" ;
    my $cigarH = &SeqAlnSunhh::parseCigar( $ta[5] );
    my $global_mismat = 0;
    for my $tk (qw/Slen Hlen/) { defined $cigarH->{$tk} and $global_mismat += $cigarH->{$tk}; }
    for my $tk (@ta[11..$#ta]) { $tk =~ m!^NM:i:(\d+)$! or next; $global_mismat += $1; last; }
    if ($prevID ne $rdID) {
      if ($prevID ne '') {
        print STDOUT "$prevID\t$prevLen\t$prevMis\n";
      }
      $prevID = $rdID;
      $prevLen = $cigarH->{'RdLen'};
      $prevMis = $global_mismat;
    } else {
      $prevMis > $global_mismat and $prevMis = $global_mismat;
    }
  }
  $prevID eq '' or print STDOUT "$prevID\t$prevLen\t$prevMis\n";
  $prevID = '';
} else {
  my (@rdStorage, %rdID2idx);
  while (<>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[5] eq '*' and next;
    defined $flag_UN{$ta[1]} and next; # Read is unmapped.
    my $rdID = (defined $flag_P{$ta[1]}) ? "$ta[0]\tR1" : "$ta[0]\tR2" ;
    my $cigarH = &SeqAlnSunhh::parseCigar( $ta[5] );
    defined $rdID2idx{$rdID} or do { push(@rdStorage, [$rdID, $cigarH->{'RdLen'}, $cigarH->{'RdLen'}]); $rdID2idx{$rdID} = $#rdStorage; };
    my $global_mismat = 0;
    for my $tk (qw/Slen Hlen/) { defined $cigarH->{$tk} and $global_mismat += $cigarH->{$tk}; }
    for my $tk (@ta[11..$#ta]) { $tk =~ m!^NM:i:(\d+)$! or next; $global_mismat += $1; last; }
    $rdStorage[ $rdID2idx{$rdID} ][2] > $global_mismat and $rdStorage[ $rdID2idx{$rdID} ][2] = $global_mismat;
  }
  for (@rdStorage) {
    print STDOUT join("\t", @$_)."\n";
  }
}

