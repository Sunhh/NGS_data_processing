#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 final.wiRepre.fmt.PAV.1.rmOlap final.wiRepre.fmt.PAV > final.wiRepre.fmt.PAV.rmOlap\n";

my $f1_olap = shift;
my $f2_pav  = shift;

my $n1 = 3;   # column-N for gene number
my $n2 = 7;   # column-N for gene IDs.
my $n3 = 11;  # column-N for PAV;

my %gene2grp = &load_rmolap($f1_olap); # {'grpID'}{$grpID}=1, {'gene2grpID'}[specIdx]{$geneID} = $grpID; {'specN'} = species_number;

my $fh2 = &openFH($f2_pav, '<');
my (@o1, %grp2idx);
while (<$fh2>) {
  chomp;
  my @ta=&splitL("\t", $_);
  if ($ta[0] eq 'groupID') {
    print STDOUT "$_\n";
    next;
  }
  if (defined $grp2idx{$ta[0]}) {
    &update_o1(\@o1, \@ta, $grp2idx{$ta[0]});
  } elsif (defined $gene2grp{'grpID'}{$ta[0]}) {
    push(@o1, [@ta]);
    $grp2idx{$ta[0]} = $#o1;
  } else {
    my $newGrpID = '';
    for (my $i=1;$i<=$gene2grp{'specN'};$i++) {
      for my $gid (split(/;/,$ta[$i+$n2-1])) {
        $gid =~ s!\s!!g;
        $gid =~ m!^\s*$! and next;
        defined $gene2grp{'gene2grpID'}[$i]{$gid} or &stopErr("[Err] Failed to find groupID for gene [$gid]\n");
        $newGrpID = $gene2grp{'gene2grpID'}[$i]{$gid};
        last;
      }
      $newGrpID ne '' and last;
    }
    $newGrpID ne '' or &stopErr("[Err] Failed to find groupID for line: $_\n");
    if (defined $grp2idx{$newGrpID}) {
      my $idx = $grp2idx{$newGrpID};
      &update_o1(\@o1, \@ta, $idx);
    } else {
      push(@o1, [$newGrpID, @ta[1..$#ta]]);
      $grp2idx{$newGrpID} = $#o1;
    }
  }
}
close($fh2);
for (@o1) {
  print STDOUT join("\t", @$_)."\n";
}

sub update_o1 {
  my ($o1H, $taH, $idx) = @_;
  defined $gene2grp{'grpID'}{$taH->[0]} and $o1H->[$idx][1] = $taH->[1];
  $o1H->[$idx][2] += $taH->[2];
  for (my $j=1;$j<=$gene2grp{'specN'};$j++) {
    $o1H->[$idx][$j+$n1-1] += $taH->[$j+$n1-1];
    $o1H->[$idx][$j+$n2-1] = join(";", grep {$_ ne ""} ( $taH->[$j+$n2-1], $o1H->[$idx][$j+$n2-1] ) );
  }
  for (my $j=$n3;$j<@$taH;$j++) {
    $taH->[$j] eq 'P' and $o1H->[$idx][$j] = 'P';
  }
  return;
}# update_o1()

sub load_rmolap {
  my ($fn) = @_;
  my $fh = &openFH($fn, '<');
  my @back;
  my %finalGrpID;
  my %b1;
  while (<$fh>) {
    chomp;
    my @ta=split(/\t/, $_);
    $finalGrpID{$ta[0]} = 1;
    for (my $i=1; $i<@ta; $i++) {
      for my $gid (split(/;/,$ta[$i])) {
        $back[$i] //= {};
        $gid =~ s!\s!!g;
        $back[$i]{$gid} = $ta[0];
      }
    }
    $b1{'specN'} //= $#ta;
    $b1{'specN'} < $#ta and $b1{'specN'} = $#ta;
  }
  close($fh);
  $b1{'grpID'} = \%finalGrpID;
  $b1{'gene2grpID'} = \@back;
  return(%b1);
}# load_rmolap()
