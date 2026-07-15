#!/bin/bash
# Runnable WGCNA demo on the bundled input/ data. To finish quickly it subsets to the top
# ~3000 most-variable genes; drop the subset step to run on the full matrix (a few minutes).
# Needs Rscript with WGCNA (+ pheatmap, RColorBrewer). Outputs go to ./_out/.
set -e
cd "$(dirname "$0")"
O=_out; mkdir -p "$O"

# input: gunzip the expression matrix, then keep the top-variable genes for speed
zcat input/dat1_rpkmMean.gz > "$O/expr.full"
perl -e 'open F,"<",$ARGV[0]; $h=<F>; my @r;
  while(<F>){chomp; my @a=split/\t/; my $id=shift @a; my($s,$ss,$n)=(0,0,0);
    for(@a){next if $_ eq "" || $_ eq "-"; $s+=$_;$ss+=$_*$_;$n++} next if $n<2;
    push @r,[($ss-$s*$s/$n)/($n-1),$_];}
  @r=sort{$b->[0]<=>$a->[0]}@r; open O,">",$ARGV[1]; print O $h;
  for my $i (0..2999){last if $i>$#r; print O $r[$i][1],"\n"} ' "$O/expr.full" "$O/expr"

# (2) build the network (defaults: log2 + bicor + signed)
Rscript run_wgcna.r --expr "$O/expr" --outdir "$O/wgcna" --prefix dat1

# (3) associate modules with the (binary) traits
Rscript down_phenoAssoc.r --network "$O/wgcna/dat1_network_ForPhenoAssoc.RData" \
  --pheno input/dat1_pheno --out "$O/wgcna/dat1" --corType bicor --binaryTrait

echo "Done. Key outputs: $O/wgcna/dat1_KME.txt  $O/wgcna/dat1.module_trait.tsv (+ pdfs)"
