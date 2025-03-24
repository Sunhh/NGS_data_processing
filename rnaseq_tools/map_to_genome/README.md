# Pipeline to do RNA-seq analysis by mapping clean reads to a single reference genome.

## (1) Map reads to a reference genome with two-step HiSat2 (No GFF files used).
- Build HISAT2 database.
```sh
hisat2-build in_chr.fa ref_db
```

- Map reads to HISAT2 database with two-step mapping.
 - Prepare an input read list with the format as the example below. dataPrefix and Fq names must be unique. Fq files can be gz.

| #sample | readGroup | library | dataPrefix | inFq1 | inFq2 | PL | PU | Others |
|---------|-----------|---------|------------|-------|-------|----|----|--------|
| S1      | S1\_Rep1  |S1\_Rep1 | S1\_Rep1   | path1 | NA    | NA | NA | NA     |
| S1      | S1\_Rep2  |S1\_Rep2 | S1\_Rep2   | path1 | NA    | NA | NA | NA     |


