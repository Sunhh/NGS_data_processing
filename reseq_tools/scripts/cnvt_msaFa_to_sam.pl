#!/usr/bin/perl
use strict;
use warnings;
use fastaSunhh;

!@ARGV and die "perl $0 in_msa_aln.fa  ref_ID > in_msa_aln.sam\n -- Convert MSA.fasta file to SAM file for IGV.\n";

my $inFa  = shift;
my $refID = shift;

my $fs=fastaSunhh->new();
my %seqs = %{$fs->save_seq_to_hash('faFile' => $inFa)};
my (@refArr, @rdIDs, %rdArrH);
for (sort {$seqs{$a}{'Order'} <=> $seqs{$b}{'Order'} } keys %seqs) {
  $refID //= $_;
  $seqs{$_}{'seq'} =~ s!\s!!g;
  if ($_ eq $refID) {
    @refArr = split(//, $seqs{$_}{'seq'});
  } else {
    push(@rdIDs, $_);
    $rdArrH{$_} = [split(//, $seqs{$_}{'seq'})];
  }
}

# Label refArr.
my @refArrPos;
{
my $prevP=0;
for (my $i=0; $i<@refArr; $i++) {
  if ($refArr[$i] eq '-') {
    $refArrPos[$i] = $prevP;
  } else {
    $prevP ++; 
    $refArrPos[$i] = $prevP;
  }
}
print STDOUT "\@HD\tVN:1.6\tSO:coordinate\n";
print STDOUT "\@SQ\tSN:$refID\tLN:$prevP\n";
print STDOUT "\@RG\tID:MSA\n";
print STDOUT "\@PG\tID:MSA_to_SAM\n";
}

# Generate rd cigar and start position.
my %rdCigar;
for my $rdID (@rdIDs) {
  my @rdArr = @{$rdArrH{$rdID}};
  my @cigar = ([0, 'M']);
  my $rdTrim = 0;
  my $mapStart = 0; # 1-based
  for (my $i=0; $i<@refArr; $i++) {
    if ($refArr[$i] eq '-') {
      # Gap in REF.
      if ($mapStart == 0) {
        # Left rd not on REF.
        $rdArr[$i] ne '-' and $rdTrim ++;
        next;
      }
      # Record insertions.
      if ($rdArr[$i] ne '-') {
        if ($cigar[-1][1] eq 'I') {
          $cigar[-1][0] ++;
        } else {
          push(@cigar, [1, 'I']);
        }
      }
    } elsif ($rdArr[$i] eq '-') {
      # Base in REF but 'D' in rd.
      $mapStart == 0 and next;
      # Record deletions.
      if ($cigar[-1][1] eq 'D') {
        $cigar[-1][0] ++; # Extend 'D'
      } else {
        push(@cigar, [1, 'D']); # Add new 'D'.
      }
    } else {
      # Base in both REF and rd. No matter if they are the same.
      if ($mapStart == 0) {
        $mapStart = $refArrPos[$i]; # First base-to-base position.
      }
      # Record matches.
      if ($cigar[-1][1] eq 'M') {
        $cigar[-1][0] ++; # Extend 'M'
      } else {
        push(@cigar, [1, 'M']); # Add new 'M'
      }
    }
  }
  # If the filling element is not used, there must be something wrong!
  $cigar[0][0] == 0 and $cigar[0][1] eq 'M' and die "Error1: rdTrim=$rdTrim;mapStart=$mapStart;firstCigar=$cigar[0][0]$cigar[0][1];\n";
  $mapStart > 0 or die "Error2: mapStart=$mapStart;\n";
  # Generate CIGAR with soft clips.
  $cigar[-1][1] eq 'I' and $cigar[-1][1] = 'S';
  $rdTrim > 0 and unshift(@cigar, [$rdTrim, 'S']);
  my $cigar_text = join("", map {"$_->[0]$_->[1]"} @cigar);
  # Generate read sequence.
  my $rdSeq = join("", grep {$_ ne "-"} @rdArr);
  ### Check consistency
  my $sum1 = 0; for (map {$_->[0]} grep {$_->[1] ne 'D'} @cigar) { $sum1 += $_; }
  my $sum2 = length($rdSeq);
  $sum1 == $sum2 or die "Error3: sum1=$sum1;sum2=$sum2\n";
  # Output SAM lines.
  print STDOUT join("\t", $rdID, 0, $refID, $mapStart, 255, $cigar_text, "*", 0, 0, $rdSeq, "*")."\n";
}

