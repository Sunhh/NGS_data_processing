#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "upDist:i", # 2000
  "downDist:i", # 0
  "refGff:s", # required
  "qryGff:s", # required
);

$opts{'upDist'}   //= 2000;
$opts{'downDist'} //= 0;

my $htxt = <<HH;
################################################################################
perl $0 -refGff prot.chr.gff3 -qryGff  struct.gff3 > struct_related_gene.list

 -upDist      [2000]
 -downDist    [0]
################################################################################
HH

for my $k1 (qw/refGff qryGff/) {
  defined $opts{$k1} or &LogInforSunhh::usage($htxt);
}

my (%locGene, %locCDS, %locUp, %locDown);
my %cdsSE;

{
  # Process reference gff3 file.
  my $rGfh = &openFH($opts{'refGff'}, '<');
  &tsmsg("[Msg] Loading refGff [$opts{'refGff'}]\n");
  while (<$rGfh>) {
    m!^\s*(#|$)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    if ($ta[2] =~ m!^CDS$!i) {
      $ta[8] =~ m!(?:^|\s|;)Parent=([^\s;]+)!i or &stopErr("[Err] Bad 1: $_\n");
      my $pid = $1;
      $cdsSE{$pid} //= [@ta[3,4,6,0]];
      scalar(@{$cdsSE{$pid}}) == 0 and $cdsSE{$pid} = [@ta[3,4,6,0]];
      $cdsSE{$pid}[0] > $ta[3] and $cdsSE{$pid}[0] = $ta[3];
      $cdsSE{$pid}[1] < $ta[4] and $cdsSE{$pid}[0] = $ta[4];
      $cdsSE{$pid}[3] eq $ta[0] or &stopErr("[Err] Bad 3: $cdsSE{$pid}[3] eq $ta[0] : $_\n");
      $locCDS{$ta[0]} //= [];
      for my $p ($ta[3]..$ta[4]) {
        $locCDS{$ta[0]}[$p]{$pid} = 1;
      }
    } elsif ( $ta[2] =~ m!^(mRNA|transcript)$!i ) {
      $ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)!i or &stopErr("[Err] Bad 2: $_\n");
      my $id = $1;
      $cdsSE{$id} //= [];
    }
  }
  close($rGfh);
  &tsmsg("[Msg] Filling the variables.\n");
  for my $mid (keys %cdsSE) {
    scalar(@{$cdsSE{$mid}}) == 0 and do { &tsmsg("[Wrn] No CDS found for mRNA [$mid]\n"); next;};
    $locGene{$cdsSE{$mid}[3]} //= [];
    for my $p ($cdsSE{$mid}[0] .. $cdsSE{$mid}[1]) { $locGene{$cdsSE{$mid}[3]}[$p]{$mid} = 1; }
    $locUp{$cdsSE{$mid}[3]} //= [];
    $locDown{$cdsSE{$mid}[3]} //= [];
    my ($u_s, $u_e, $d_s, $d_e);
    if ($cdsSE{$mid}[2] eq '-') {
      $u_s = $cdsSE{$mid}[1]+1;
      $u_e = $cdsSE{$mid}[1]+$opts{'upDist'};
      $d_s = $cdsSE{$mid}[0]-$opts{'downDist'};
      $d_e = $cdsSE{$mid}[0]-1;
      $d_s > 0 or $d_s = 1;
      for my $p ($u_s .. $u_e) { $locUp{$cdsSE{$mid}[3]}[$p]{$mid} //= 1; }
      for my $p ($d_s .. $d_e) { $locDown{$cdsSE{$mid}[3]}[$p]{$mid} //= 1; }
    } elsif ($cdsSE{$mid}[2] eq '+' or $cdsSE{$mid}[2] eq '.') {
      $u_s = $cdsSE{$mid}[0]-$opts{'upDist'};
      $u_e = $cdsSE{$mid}[0]-1;
      $u_s > 0 or $u_s = 1;
      $d_s = $cdsSE{$mid}[1]+1;
      $d_e = $cdsSE{$mid}[1]+$opts{'downDist'};
      for my $p ($u_s .. $u_e) { $locUp{$cdsSE{$mid}[3]}[$p]{$mid} //= 1; }
      for my $p ($d_s .. $d_e) { $locDown{$cdsSE{$mid}[3]}[$p]{$mid} //= 1; }
    } else {
      &stopErr("[Err] Bad 4: @{$cdsSE{$mid}}\n");
    }
  }
}

# Return: ID/Name/var_len/type2
sub parse_sv {
  my ($str) = @_;
  my %back;
  $str =~ m!ID=([^\s;]+);Name=([^\s;]+)! or die "bad ID: $str\n";
  $back{'ID'}   = $1;
  $back{'Name'} = $2;
  $str =~ m!;(overlap|subst|del|ins|blk)_len=(\d+)! and $back{'var_len'} = $2;
  $back{'var_len'} //= -1;
  if ($back{'Name'} eq 'substitution' and $back{'var_len'} == 1) {
    $back{'type2'} = 'SNP';
  } elsif ($back{'Name'} =~ m!gap|ATGCN!) {
    $back{'type2'} = 'GAP';
  } elsif ($back{'Name'} =~ m!^(deletion|collapsed_repeat|collapsed_tandem_repeat)$!) {
    $back{'type2'} = 'DELETION';
  } elsif ($back{'Name'} =~ m!^(insertion|duplication|tandem_duplication|unaligned_end|unaligned_beginning|relocation-insertion|translocation\-insertion)$!) {
    $back{'type2'} = 'INSERTION';
  } elsif ($back{'Name'} =~ m!^inversion$!) {
    $back{'type2'} = 'INVERSION';
  } elsif ($back{'Name'} =~ m!^reshuffling\-part_!) {
    $back{'type2'} = "RESHUFFLING";
  } elsif ($back{'Name'} =~ m!^(translocation|relocation)!) {
    $back{'type2'} = uc($1);
  } else {
    $back{'type2'} = uc($back{'Name'});
  }
  return(%back);
}# parse_sv ()
{
  # Find affected genes.
  my $qGfh = &openFH($opts{'qryGff'}, '<');
  &tsmsg("[Msg] Processing qryGff [$opts{'qryGff'}]\n");
  print STDOUT join("\t", qw/VAR_ID  Chr_ID Chr_start Chr_end Affected_type Affected_genes VAR_type1 VAR_type2 VAR_len VAR_annot/)."\n";
  while (<$qGfh>) {
    m!^\s*(#|$)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    my %info = &parse_sv($ta[8]);
    # $ta[8] =~ m!(?:^|;|\s)ID=([^\s;]+)!i or &stopErr("[Err] Bad 5:$_\n");
    my $id = $info{'ID'};
    my ($chr, $s, $e) = @ta[0,3,4];
    my $hasFound = 0;
    my %effGenes;
    # Check CDS
    for my $p ($s .. $e) {
      defined $locCDS{$chr}[$p] or next;
      $hasFound = 1;
      for my $mid (keys %{$locCDS{$chr}[$p]}) {
        $effGenes{$mid}++;
      }
    }
    if ($hasFound == 1) {
      print STDOUT join("\t", $id, $chr, $s, $e, 'CDS', join(";", sort { $effGenes{$b} <=> $effGenes{$a} } keys %effGenes), $info{'type2'}, $info{'Name'}, $info{'var_len'}, $ta[8])."\n";
      next;
    }
    # Check gene
    for my $p ($s .. $e) {
      defined $locGene{$chr}[$p] or next;
      $hasFound = 1;
      for my $mid (keys %{$locGene{$chr}[$p]}) {
        $effGenes{$mid}++;
      }
    }
    if ($hasFound == 1) {
      print STDOUT join("\t", $id, $chr, $s, $e, 'intron', join(";", sort { $effGenes{$b} <=> $effGenes{$a} } keys %effGenes), $info{'type2'}, $info{'Name'}, $info{'var_len'}, $ta[8])."\n";
      next;
    }
    # Check Up
    for my $p ($s .. $e) {
      defined $locUp{$chr}[$p] or next;
      $hasFound = 1;
      for my $mid (keys %{$locUp{$chr}[$p]}) {
        $effGenes{$mid}++;
      }
    }
    ### Select Up
    if ($hasFound == 1) {
      my @gidsP = sort { $cdsSE{$a}[0] <=> $cdsSE{$b}[0] } grep { $cdsSE{$_}[2] eq '+' } keys %effGenes;
      my @gidsM = sort { $cdsSE{$b}[1] <=> $cdsSE{$a}[1] } grep { $cdsSE{$_}[2] eq '-' } keys %effGenes;
      my @gidsFinal;
      scalar(@gidsP) > 0 and push(@gidsFinal, $gidsP[0]);
      scalar(@gidsM) > 0 and push(@gidsFinal, $gidsM[0]);
      print STDOUT join("\t", $id, $chr, $s, $e, 'upstream', join(";", @gidsFinal), $info{'type2'}, $info{'Name'}, $info{'var_len'}, $ta[8])."\n";
      next;
    }
    # Check Down
    for my $p ($s .. $e) {
      defined $locDown{$chr}[$p] or next;
      $hasFound = 1;
      for my $mid (keys %{$locDown{$chr}[$p]}) {
        $effGenes{$mid}++;
      }
    }
    ### Select Down
    if ($hasFound == 1) {
      my @gidsP = sort { $cdsSE{$b}[1] <=> $cdsSE{$a}[1] } map { $cdsSE{$_}[2] eq '+' } keys %effGenes;
      my @gidsM = sort { $cdsSE{$a}[0] <=> $cdsSE{$b}[0] } map { $cdsSE{$_}[2] eq '-' } keys %effGenes;
      my @gidsFinal;
      scalar(@gidsP) > 0 and push(@gidsFinal, $gidsP[0]);
      scalar(@gidsM) > 0 and push(@gidsFinal, $gidsM[0]);
      print STDOUT join("\t", $id, $chr, $s, $e, 'downstream', join(";", @gidsFinal), $info{'type2'}, $info{'Name'}, $info{'var_len'}, $ta[8])."\n";
      next;
    }
    # Rest.
    print STDOUT join("\t", $id, $chr, $s, $e, "-", "-", $info{'type2'}, $info{'Name'}, $info{'var_len'}, $ta[8])."\n";
  }
  close($qGfh);
}

