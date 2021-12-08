#!/usr/bin/perl
# 20211207: Currently I assume all the sequences are included with specific naming rules.
#   The first column is sequence ID, the second column is the index used in later assembly, and the third is the length of current segment. 
#   For example: 
#     (1) The three lines below mean the total length of ptg000011l is 2697931 (1382931+5000+1310000), and the second sequence (ptg000011l:::fragment_2:::debris) is dropped as debris. 
#     >ptg000011l:::fragment_1 14 1382931
#     >ptg000011l:::fragment_2:::debris 15 5000
#     >ptg000011l:::fragment_3 16 1310000
#     (2) The two lines below tell that the total length of ptg000101l is 133473 (0+133473), and the first sequence sizes zero bp. I don't know why it occurs, but it doesn't matter much. 
#     >ptg000101l:::fragment_1:::debris 130 0
#     >ptg000101l:::fragment_2 131 133473
#     (3) ">hic_gap_375 379 500" means this is an additional gap sequence with 500 Ns. 
#   There is something wrong in the revised .assembly file too, for example: 
#     (1) The gap is incorrectly left and failed to inserted to the new scaffold joint, so I'd like to ignore the gap information and add gaps on my own. 
#     # C31_FINAL (before manual revision)
#     >hic_gap_375 375 500
#     >ptg000065l 16 7690470
#     >ptg000002l 11 9987887
#     -16 375 11
#     # C31_F1 (after manual revision)
#     >ptg000065l 20 7690470
#     >ptg000002l 11 9987887
#     >hic_gap_375 379 500
#     20 -11 379
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and &stopErr("perl $0 C31_F1.assembly > C31_F1.agp\n"); 

my $infer_gap_len = 1; 
my $gap_len = 500; 
my $asm_prefix = "HiC_scaffold_"; 

my (@segments); 
my %toResolve; 
my %ngapIdx; 
my @assemblies; 
while (<>) {
  chomp; 
  if (m!^\>(\S+)\s+(\d+)\s+(\d+)$!) {
    my ($name0, $i0, $l0) = ($1, $2, $3); 
    $name0 =~ s!:::debris$!!; 
    $segments[$i0] = [$name0, 1, $l0]; 
    if      ($name0 =~ m!^(\S+):::fragment_(\d+)$!) {
      my ($name1, $idx1) = ($1, $2); 
      push(@{$toResolve{$name1}}, [$idx1, $name0, $i0, $l0]); 
    } elsif ($name0 =~ m!^hic_gap_\d+$!) {
      $ngapIdx{$i0} = 1; 
    }
    
  } elsif (m!^(\-?\d+)(?:\s+\-?\d+)*$!) {
    my @ta = split(/\s+/, $_); 
    my @tasm; 
    for my $tb (@ta) {
      my ($ti, $ts); 
      if ($tb =~ m!^\-(\d+)$!) {
        $ti = $1; $ts = "-"; 
      } elsif ($tb =~ m!^(\d+)$!) {
        $ti = $1; $ts = "+"; 
      } else {
        &stopErr("[Err] Bad assembly index in line: $_\n"); 
      }
      defined $ngapIdx{$ti} and next; 
      push(@tasm, [$ti, $ts]); 
    }
    push(@assemblies, [@tasm]); 
  } else {
     &stopErr("[Err] Bad line: $_\n"); 
  }
}

# Resolve fragments's positions. 
for my $name1 (keys %toResolve) {
  my @ta = sort { $a->[0] <=> $b->[0] } @{$toResolve{$name1}}; 
  my $ttl_len = 0; 
  for my $tb (@ta) {
    $segments[$tb->[2]] = [ $name1, $ttl_len+1, $ttl_len+$tb->[3] ]; 
    $ttl_len += $tb->[3]; 
  }
}

# Resolve the gap length
if ($infer_gap_len == 1 and scalar(keys %ngapIdx) == 1) {
  my ($i0) = keys %ngapIdx; 
  $gap_len = $segments[$i0][2]; 
}

# Build the AGP file. 
### Resolve the structure and length of assemblies. 
my @asm; 
for my $ta (@assemblies) {
  push(@asm, [0, []]); # [assembly_length, assembly_structure]
  my $agp_idx = 0; 
  for (my $i=0; $i<@$ta; $i++) {
    my ($i0, $str) = @{$ta->[$i]}; 
    my $seg_len = $segments[$i0][2]-$segments[$i0][1]+1; 
    $agp_idx ++; 
    push(@{$asm[-1][1]}, [$asm[-1][0]+1, $asm[-1][0]+$seg_len, $agp_idx, "W", $segments[$i0][0], $segments[$i0][1], $segments[$i0][2], $str]); 
    $asm[-1][0] += $seg_len; 
    if ($i < $#$ta) {
      $agp_idx ++; 
      push(@{$asm[-1][1]}, [$asm[-1][0]+1, $asm[-1][0]+$gap_len, $agp_idx, "U", $gap_len, "scaffold", "yes", "paired-ends"]); 
      $asm[-1][0] += $gap_len; 
    }
  }
}# End for my $ta (@assemblies) 
### Add names. 
@asm = sort { $b->[0] <=> $a->[0] } @asm; 
for (my $i=0; $i<@asm; $i++) {
  my $j = $i + 1; 
  my $out_seqname = ${asm_prefix}.$j; 
  for my $ta (@{$asm[$i][1]}) {
    print join("\t", $out_seqname, @$ta)."\n"; 
  }
}
