# Scripts to process resequencing data.

## Workflow to count read pairs spanning target region.
It is recommended to use whole genome as the reference for alignment, because truncated reference can force repeat-related read pairs aligned to the junction site.
- This method can be used to genotype large deletions/insertions by appropriate target locatioin.
  - Tandem duplication (insertion): targeting junction site between two duplicates.
  - Large deletion: targeting read pairs spanning the entire deleted region.
  - Large exogenous insertion: targeting read pairs spanning the two ends of inserted region.

```sh
# BWA-MEM alignments. "bwa mem -M | samtools view -h -F 12| samtools sort -o PE2allele.bam"
# It is recommended to extend junction site 10-bp outside to avoid false mapping at the read ends forced by aligners.
# Nearby S and E are used as positive control to test if read depth is enought to capture junction site.
perl scripts/cntRd_spanJunctionSite_inBam.pl -in_bam PE2allele.bam -target_loc junction_seq:junction_S-junction_E  -positive_control_loc junction_seq:nearby_E-nearby_E > PE2allele.read_counts
```

There are three columns in file PE2allele.read\_counts: (0) input file name, (1) number of pairs spanning target, (2) number of pairs spanning control.

## Workflow to genotype large deletions in aligned BAM file.
- Input files:
  - List of long deletion variants: Example file `list.long_deletions` with columns of `variant ID`, `chromosome`, `deletion_start`, `deletion_end`;
  - List of accession to bam file mapping: Example file `list.sample_bam` with columns of `sample ID`, `bam file path`;
  - BAM alignments: Example file like `ARO19494.dedup.bam`;
- Result files:
  - `out_geno-melt.tab`: `variant ID`, `Sample ID`, `0/0 or 1/1` genotype, read pairs spanning left boundary, read pairs spanning right boundary.
  - `out_geno-mat.tab` : Line 1: `Sample`, `variant IDs`; Line 2 and following: `Sample ID`, `0/0 or 1/1` genotypes;

```sh
perl genotype_longDEL_byListBam.pl  list.long_deletions  list.sample_bam  out_geno
perl cnvt_melt_to_matrix.pl out_geno-melt.tab > out_geno-mat.tab
```

## File format convertioin.
- Convert ClustalW output MSA.fasta file into a SAM file.

```sh
perl cnvt_msaFa_to_sam.pl  in_MSA.fasta reference_ID_in_MSA  > out_MSA.sam
```
