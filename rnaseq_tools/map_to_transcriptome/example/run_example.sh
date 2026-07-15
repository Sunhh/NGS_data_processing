#!/bin/bash
# Demo of the transcriptome DEG steps (3->5) on the toy kallisto data in this folder.
# Step 2 (salmon/kallisto quant) needs the binaries + real reads, so this starts from
# per-sample quant dirs (kallisto_out/*/abundance.tsv). Salmon is identical: point at
# salmon_out/*/ and pass -t salmon. Needs: Rscript (tximport, DESeq2), perl.
set -e
cd "$(dirname "$0")"
P=..; O=_out; mkdir -p "$O"

# (3) merge per-sample kallisto quant -> gene(OG)-level tables + RDS
ls -d "$PWD"/kallisto_out/*/ | sed 's#/$##' > "$O/list.sample_dir"
Rscript $P/merge_tx_quant_to_gene.r -s "$O/list.sample_dir" -g tx2gene.tsv -t kallisto --quant_file abundance.tsv -o "$O/comb"

# (4) tximport DESeq2 -> FDR-only table (ds.<g1>_VS_<g2> per comparison)
Rscript $P/run_deseq2_tximport.r -r "$O/comb.RDS" -m list.sample_meta -c list.comparison -o "$O/joint-fdr.txt"

# (5) DEG calls in vClaude (shared caller; same list.comparison as --pair)
perl $P/../DEG_byList_vClaude.pl --tpm "$O/comb.gene_tpm.tsv" --pair list.comparison --fdr "$O/joint-fdr.txt" \
     --fdr_cut 0.05 --fc 2 --out "$O/DEG_calls.tsv"

echo "Done. Outputs in $O/ : comb.gene_tpm.tsv  joint-fdr.txt  DEG_calls.tsv"
