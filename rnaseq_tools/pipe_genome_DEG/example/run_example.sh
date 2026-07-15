#!/bin/bash
# Runnable demo of the DEG steps (2->5) on the toy inputs in this folder.
# Steps 0-1 (rRNA removal, HISAT2, featureCounts) need real reads + external tools,
# so this script starts from the featureCounts-style count matrix joint-cnt_rename.txt.
# Needs: perl, Rscript with DESeq2 + edgeR. Outputs go to ./_out/.
set -e
cd "$(dirname "$0")"
P=..                       # scripts live one level up
O=_out; mkdir -p "$O"

# (2) counts -> TPM (each column sums to 1e6)
perl $P/cnvt_cnt_to_normExpr.pl --cntFn joint-cnt_rename.txt --out_type tpm > $O/joint-tpm.txt
#     ...RPKM instead:  --out_type rpkm [--libFn sample<TAB>total_mapped_reads]

# (3) FDR only, one call over all comparisons (col1 picks DESeq2 / edgeR GLM / edgeR classic)
cut -f1,7- joint-cnt_rename.txt > $O/joint-cnt_matrix.txt         # drop featureCounts annotation cols
perl $P/DEG_byList.pl -rdCntFn $O/joint-cnt_matrix.txt -compareList list.comparison -outFDRFn $O/joint-fdr.txt
#     inspect/extract the exact R instead of running it:  -writeRcode $O/Rcode

# (4) DEG calls in vClaude (adjustable --fdr_cut / --fc); reuse list.comparison as --pair
perl $P/DEG_byList_vClaude.pl --tpm $O/joint-tpm.txt --pair list.comparison --fdr $O/joint-fdr.txt \
     --fdr_cut 0.05 --fc 2 --out $O/DEG_calls.tsv

# (5) roll gene-level DEG labels up to orthogroups
perl -e 'open F,"<",$ARGV[0]; my@h=split/\t/,scalar<F>; chomp@h; my@k=(0);
  for my $i (1..$#h){push @k,$i if $h[$i]=~/\.DEG$/;} $h[0]="GeneID";
  print join("\t",@h[@k]),"\n"; while(<F>){chomp;my@a=split/\t/;print join("\t",@a[@k]),"\n";}' \
  $O/DEG_calls.tsv > $O/DEG_labels.tsv
perl $P/cnvt_gene2group_DEGlabel.pl synOG.grp $O/DEG_labels.tsv > $O/DEG_group_labels.tsv

echo "Done. Key outputs in $O/ :"
echo "  joint-tpm.txt  joint-fdr.txt  DEG_calls.tsv  DEG_group_labels.tsv"
