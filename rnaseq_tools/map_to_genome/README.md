# Pipeline to do RNA-seq analysis by mapping clean reads to a single reference genome.

## (1) Map reads to a reference genome with two-step HiSat2 (No GFF files used).
- Build HISAT2 database.
```sh
hisat2-build in_chr.fa ref_db
```

- Map reads to HISAT2 database with two-step mapping.
  - Prepare an input read list (`list.in_rd`) with the format as the example below. dataPrefix and Fq names must be unique. Fq files can be gz.

| #sample | readGroup | library | dataPrefix | inFq1 | inFq2 | PL | PU | Others |
|---------|-----------|---------|------------|-------|-------|----|----|--------|
| S1      | S1\_Rep1  |S1\_Rep1 | S1\_Rep1   | path1 | NA    | NA | NA | NA     |
| S1      | S1\_Rep2  |S1\_Rep2 | S1\_Rep2   | path1 | path2 | NA | NA | NA     |

 - Run 2-pass HISAT2 alignments.

```sh
mkdir out_bams/
perl runHisat2_with2pass.pl \
 -cpuN         20 \
 -in_pref_list list.in_rd \
 -db_hisat2    ref_db \
 -wrk_dir      out_bams/ \
 -para_hisat2  ' -p 4 --dta --dta-cufflinks -q --phred33 --rna-strandness R '

```


