#!/usr/bin/perl
# 11/17/2021 : Edit for options. 
# 04/01/2022 : Fix the bug when there are two exactly the same contigs in the file, both of them will be removed.
#   Use more codes instead of invoking other scripts.
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
  "min_ident:f", # 0-1
  "min_cover:f", # 0-1
  "maxLen2Chk:i", # 300e3
  "bn_eval:f",    # 1e-10
  "bn_wordsize:i", # 100
  "bn_ident:f",    # 0
  "bn_wordsizeL2L:i", # 1000
  "skip_long2long!",
  "help!"
); 

$opts{'min_ident'} //= 0.99; 
$opts{'min_cover'} //= 0.99; 
$opts{'maxLen2Chk'}     //= 300e3; 
$opts{'bn_wordsize'}    //= 100; 
$opts{'bn_wordsizeL2L'} //= 1e3;
$opts{'bn_eval'}     //= 1e-10; 
$opts{'bn_ident'}    //= 0; 

my $min_identR = $opts{'min_ident'}; 
my $min_covR   = $opts{'min_cover'}; 

my $max_seqlen  = $opts{'maxLen2Chk'}; # 300kb is a good length maximum to remove redundancy. 
my $bn_wordsize = $opts{'bn_wordsize'};
my $bn_wordsizeL2L = $opts{'bn_wordsizeL2L'};
my $bn_eval     = $opts{'bn_eval'}; 
my $bn_ident    = $opts{'bn_ident'}; 

my $help_txt = <<HH; 
perl $0 X20.ctg.fa

-min_ident   [$min_identR]
-min_cover   [$min_covR]
-maxLen2Chk  [$max_seqlen]
-bn_wordsize [$bn_wordsize]
-bn_eval     [$bn_eval]
-bn_ident    [$bn_ident]
-bn_wordsizeL2L [$bn_wordsizeL2L]
-skip_long2long 

HH


!@ARGV and die $help_txt; 

my $fn = shift; 
my @torm; 

# Step 1: Remove short contigs (<= $max_seqlen) that are covered by any longer contigs.
my %seqs = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fn ) };
my %grp;
open O1,'>',"$fn.short.fa" or &stopErr("[Err] open [$fn.short.fa]\n");
for (keys %seqs) {
  $seqs{$_}{'seq'} =~ s!\s!!g;
  $seqs{$_}{'len'} = length($seqs{$_}{'seq'});
  if ($seqs{$_}{'len'} >= $max_seqlen) {
    $grp{'long'}{$_} = $seqs{$_}{'len'};
  } else {
    $grp{'short'}{$_} = $seqs{$_}{'len'};
    print O1 ">$_\n$seqs{$_}{'seq'}\n";
  }
}
close O1;
push(@torm, "$fn.short.fa");


my @incl_1; # ([long_id, short_id, long_len, short_len])
{
  # Run blastn
  &runCmd("makeblastdb -dbtype nucl -in $fn"); # I don't want to remove short sequences in the database because there might be short primary contigs that covers shorter redundancy. 
  for my $suff (qw/not ntf nto ndb nhr nin nsq/) {
    push(@torm, "$fn.$suff");
  }
  my $cmd = "blastn -task megablast ";
  $cmd .= " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand' ";
  $cmd .= " -num_threads 50 -word_size $bn_wordsize -perc_identity $bn_ident -evalue $bn_eval ";
  $cmd .= " -max_hsps 5 -max_target_seqs 50 -dust no ";
  $cmd .= " -db $fn ";
  $cmd .= " -query $fn.short.fa ";
  $cmd .= " -out $fn.short2All.bn6 ";
  &runCmd($cmd);
  push(@torm, "$fn.short2All.bn6");
  # Get relationship
  open F2,'<',"$fn.short2All.bn6" or &stopErr("[Err] open [$fn.short2All.bn6]\n");
  while (<F2>) {
    chomp;
    my @ta = split(/\t/, $_);
    $ta[0] eq $ta[1] and next;
    $ta[2] >= 100 * $min_identR or next;
    ($ta[7]-$ta[6]+1) >= $min_covR * $ta[12] or next;
    if      ($ta[12] > $ta[13]) {
      push(@incl_1, [@ta[0,1,12,13]]);
    } elsif ($ta[12] < $ta[13]) {
      push(@incl_1, [@ta[1,0,13,12]]);
    } else {
      push(@incl_1, [@ta[0,1,12,13]]);
      push(@incl_1, [@ta[1,0,13,12]]);
    }
  }
  close F2;
  @incl_1 = sort { $a->[2] <=> $b->[2] } @incl_1;
  open O2,'>',"$fn.short2All.redID" or &stopErr("[Err] open [$fn.short2All.redID]\n");
  for (@incl_1) {
    defined $seqs{$_->[0]} or next;
    if ( defined $seqs{$_->[1]} ) {
      print O2 join("\t", $_->[1], $_->[0], $_->[3], $_->[2])."\n";
      delete $seqs{$_->[1]};
    }
  }
  close O2;
}

# Step 2: Remove long contigs (> $max_seqlen) that are covered by any longer contigs.
my @incl_2;
unless ($opts{'skip_long2long'}) {
  open O3,'>',"$fn.long.fa" or &stopErr("[Err] [$fn.long.fa]\n");
  for my $id (grep { defined $seqs{$_} } keys %{$grp{'long'}}) {
    print O3 ">$id\n$seqs{$id}{'seq'}\n";
  }
  close O3;
  push(@torm, "$fn.long.fa");
  &runCmd("makeblastdb -dbtype nucl -in $fn.long.fa");
  for my $suff (qw/not ntf nto ndb nhr nin nsq/) {
    push(@torm, "$fn.long.fa.$suff");
  }
  my $cmd = "blastn -task megablast ";
  $cmd .= " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand' ";
  $cmd .= " -num_threads 50 -word_size $bn_wordsizeL2L -perc_identity $bn_ident -evalue $bn_eval ";
  $cmd .= " -max_hsps 5 -max_target_seqs 5 -dust no ";
  $cmd .= " -db $fn.long.fa ";
  $cmd .= " -query $fn.long.fa ";
  $cmd .= " -out $fn.long2long.bn6";
  &runCmd($cmd);
  push(@torm, "$fn.long2long.bn6");
  # Get relationship
  open F3,'<',"$fn.long2long.bn6" or &stopErr("[Err] open [$fn.long2long.bn6]\n");
  while (<F3>) {
    chomp;
    my @ta = split(/\t/, $_);
    $ta[0] eq $ta[1] and next;
    $ta[2] >= 100 * $min_identR or next;
    ($ta[7]-$ta[6]+1) >= $min_covR * $ta[12] or next;
    if      ($ta[12] > $ta[13]) {
      push(@incl_2, [@ta[0,1,12,13]]);
    } elsif ($ta[12] < $ta[13]) {
      push(@incl_2, [@ta[1,0,13,12]]);
    } else {
      push(@incl_2, [@ta[0,1,12,13]]);
      push(@incl_2, [@ta[1,0,13,12]]);
    }
  }
  close F3;
  @incl_2 = sort { $a->[2] <=> $b->[2] } @incl_2;
  open O4,'>',"$fn.long2long.redID" or &stopErr("[Err] open [$fn.long2long.redID]\n");
  for (@incl_2) {
    defined $seqs{$_->[0]} or next;
    if ( defined $seqs{$_->[1]} ) {
      print O4 join("\t", $_->[1], $_->[0], $_->[3], $_->[2])."\n";
      delete $seqs{$_->[1]};
    }
  }
  close O4;
}

# Generate resulting files.
&runCmd("cat $fn.short2All.redID $fn.long2long.redID > $fn.redID");
push(@torm, "$fn.short2All.redID", "$fn.long2long.redID");
open OO,'>',"$fn.noRed.fa" or &stopErr("[Err] open [$fn.noRed.fa]\n");
for (sort { $seqs{$a}{'Order'} <=> $seqs{$b}{'Order'} } keys %seqs) {
  print OO ">$seqs{$_}{'head'}\n$seqs{$_}{'seq'}\n";
}
close OO;
&runCmd("deal_fasta.pl -N50 $fn.noRed.fa > $fn.noRed.fa.N50");

for my $a (@torm) {
	system "rm -rf $a"; 
}
