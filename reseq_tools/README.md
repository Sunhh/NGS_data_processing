# Scripts to process resequencing data.

## Workflow to count read pairs spanning target region.
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


