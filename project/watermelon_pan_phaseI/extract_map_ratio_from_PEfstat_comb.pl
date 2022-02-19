#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 comb.fstat_CLpan > comb.fstat_CLpan.tbl\n";

my @oKey =       qw/filename totalPE_rd anyMap_rdPerc totalPair_rdPerc properPair_rdPerc anyMap_rd totalPair_rd properPair_rd/;
print join("\t", qw/filename totalPE_rd anyMap_rdPerc totalPair_rdPerc properPair_rdPerc anyMap_rd totalPair_rd properPair_rd/)."\n";
for my $fn0 (@ARGV) {
  open F0,'<',$fn0 or die "[Err] $fn0\n";
  my %h;
  $h{'mode'} = 'single';
  $h{'filename'} = $fn0;
  my $cnt0 = 0;
  while (<F0>) {
    chomp;
    if (m!^:::+$!) {
      $h{'mode'} = 'multi';
      $cnt0 ++;
      if ($cnt0 % 2 == 1) {
        $cnt0 > 1 and &outRecord(\%h);
        my $fn1 = <F0>;
        chomp($fn1);
        $h{'filename'} = $fn1;
      }
      next;
    }
    if (m!^(\d+)\s.+paired in sequencing$!) {
      $h{'totalPE_rd'} = $1;
    } elsif (m!^(\d+)\s.+properly paired!) {
      $h{'properPair_rd'} = $1;
    } elsif (m!^(\d+)\s.+with itself and mate mapped$!) {
      $h{'totalPair_rd'} = $1;
    } elsif (m!^(\d+)\s.+singletons!) {
      $h{'SEmap_rd'} = $1;
    }
  }
  close F0;
  if (defined $h{'filename'}) {
    &outRecord(\%h);
  }
}

sub outRecord {
  my ($hr) = @_;
  defined $hr->{'filename'} or do { %$hr = (); return(); };
  for my $k1 (qw/totalPair_rd SEmap_rd/) {
    defined $hr->{$k1} or die "k1=$k1 in $hr->{'filename'}\n";
  }
  $hr->{'anyMap_rd'} = $hr->{'totalPair_rd'} + $hr->{'SEmap_rd'};
  $hr->{'anyMap_rdPerc'}    = sprintf("%.2f", 100 * $hr->{'anyMap_rd'} / $hr->{'totalPE_rd'} );
  $hr->{'totalPair_rdPerc'} = sprintf("%.2f", 100 * $hr->{'totalPair_rd'} / $hr->{'totalPE_rd'} );
  $hr->{'properPair_rdPerc'} = sprintf("%.2f", 100 * $hr->{'properPair_rd'} / $hr->{'totalPE_rd'} );
  for my $k1 (@oKey) {
    $hr->{$k1} //= "NA";
  }
  print STDOUT join("\t", @{$hr}{@oKey})."\n";
  %$hr = ();
  return;
}# End outRecord()

# ::::::::::::::
# BAM/ASM1001_CharlestonGray.CLpan.bam.fstat
# ::::::::::::::
# 150581930 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 1672622 + 0 supplementary
# 0 + 0 duplicates
# 150354588 + 0 mapped (99.85% : N/A)
# 148909308 + 0 paired in sequencing
# 74454654 + 0 read1
# 74454654 + 0 read2
# 136942396 + 0 properly paired (91.96% : N/A)
# 148652248 + 0 with itself and mate mapped
# 29718 + 0 singletons (0.02% : N/A)
# 7736480 + 0 with mate mapped to a different chr
# 6320922 + 0 with mate mapped to a different chr (mapQ>=5)
# ::::::::::::::
# BAM/ASM1002_97103.CLpan.bam.fstat
# ::::::::::::::
# 93147814 + 0 in total (QC-passed reads + QC-failed reads)
