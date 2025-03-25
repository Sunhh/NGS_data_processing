# Pipeline to do RNA-seq analysis by mapping clean reads to transcriptome fasta with Salmon.
- Reference website: `https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html`
- Reference guide: `https://github.com/COMBINE-lab/salmon?tab=readme-ov-file`
- Referencing videos by Ming Tang.
  - `https://www.youtube.com/watch?v=_Q8fYokTCTs&ab_channel=chatomics`
  - `https://www.youtube.com/watch?v=RWpY7EqHOUw&ab_channel=chatomics`

Currently, I am using cDNA-only index instead of SAF genome index.

- Data structure:
  - `db/`        : database and basic sequences.
  - `data/`      : information of samples and comparisons.
  - `salmon_out/`: Salmon run and output directory.
  - `result/`    : result files extracted from Salmon output that can be used further.
  - `DEGs/`      : DEG output tables.

## Group and project cDNA sequences.
- Could be orthologous groups from several genomes: `97103.c.fa`, `PI296341.c.fa`;
  - Result files
    - `slct-synOG.c.fa`: transcript sequences.
    - `slct-synOG.trans2gene.tsv`: Mapping transcripts to gene.

```sh
cat db/97103.c.fa db/PI296341.c.fa > db/slct-synOG.c.fa
deal_table.pl data/update_t2t/t2t-synOG.grp -column 0,6,46  | awk -F"\t" '!($2 == "" && $3 == "")' > data/slct-synOG.grp
perl cnvt_synOGgrp_to_trans2gene.pl data/slct-synOG.grp > data/slct-synOG.trans2gene.tsv
```

## Build Salmon index.

```sh
mamba activate salmon
salmon index -t db/slct-synOG.c.fa  -i db/synOG_index
mamba deactivate salmon
```

## Quantify gene expression using Salmon
- Quantify each sample.
  - `data/list.in_rd-fruit`: Format as used by `runHisat2_with2pass.pl`.

```sh
mkdir salmon_out/
mamba activate salmon
deal_table.pl -column 3,4 data/list.in_rd-fruit |tail -n +2|perl -e 'while (<>) {chomp;my @a=split;print "salmon quant -i db/synOG_index  -l A -p 8 -r $a[1] -o salmon_out/$a[0]\n";}' > cx1quant
nohup run_cmd_in_batch.pl -cpuN 8 -wait_sec 1 cx1quant
mamba deactivate
```
- Merge quantifications and keep a RDS file for `R/tximport`;
  - `comb.RDS`: Input file for R script `get_salmon_gene_quant_batch.r`.
  - `comb.gene_tpm.tsv`: Scaled gene-level TPM used finally.
  - `comb.gene_lengths.csv`: Averaged gene length. This value can be slightly different among samples.
  - `comb.gene_counts.tsv` : Scaled gene-level counts used in DESeq2.

```sh
mkdir result/
awk -F"\t" 'NR > 1 {print "salmon_out/"$4}' data/list.in_rd-fruit > result/list.sample_dir
Rscript get_salmon_gene_quant_batch.r  result/list.sample_dir  data/slct-synOG.trans2gene.tsv  result/comb
```

## Find DEGs using R/DESeq2.
- `list.sample_meta-fruit`: Has two columns, 'sample' and 'group'.

```sh
mkdir -p DEGs/sep/
# Find DEGs for each comparison.
perl -e 'while (<>) {chomp;my @a=split; print "Rscript run_deseq2_salmon.r -r result/comb.RDS  -m data/list.sample_meta-fruit  -b $a[0] -t $a[1] -o DEGs/sep/res-$a[2]\n";}' data/list.comparison-fruit | grep -v outPrefix > cx2deg
```

## Combine DEGs.
- `list.DEG_meta.txt`: Mapping DESeq2 result file and the comparison name.

```sh
perl map_to_genome/combine_DEGs.pl data/list.DEG_meta.txt > DEGs/combined_DEGs.txt
perl map_to_genome/label_DEGs.pl DEGs/combined_DEGs.txt > DEGs/label_DEGs.tsv
```

