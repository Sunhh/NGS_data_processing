# Pipeline to do RNA-seq analysis by mapping clean reads to a single reference genome.

## (1) Map reads to a reference genome with two-step HiSat2 (No GFF files used).
- Build HISAT2 database.
```sh
hisat2-build in_chr.fa ref_db
```

- Prepare an input read list (`list.in_rd`) with the format as the example below. dataPrefix and Fq names must be unique. Fq files can be gz.

| #sample | readGroup | library | dataPrefix | inFq1 | inFq2 | PL | PU | Others |
|---------|-----------|---------|------------|-------|-------|----|----|--------|
| S1      | S1\_Rep1  |S1\_Rep1 | S1\_Rep1   | path1 | NA    | NA | NA | NA     |
| S1      | S1\_Rep2  |S1\_Rep2 | S1\_Rep2   | path1 | path2 | NA | NA | NA     |

- Map reads to HISAT2 database with two-step mapping.
  - Resulting files: `*_fixNH.bam` files in `out_bams/`.
  - Parameter ` --rna-strandness R ` fits strand-specific RNA-seq with reads on the reverse strand.

```sh
mkdir out_bams/
perl runHisat2_with2pass.pl \
 -cpuN         20 \
 -in_pref_list list.in_rd \
 -db_hisat2    ref_db \
 -wrk_dir      out_bams/ \
 -para_hisat2  ' -p 4 --dta --dta-cufflinks -q --phred33 --rna-strandness R '

```

## (2) Count on-gene reads with featureCounts.
- Prepare GTF file from GFF3 file using gffread.
```sh
gffread  ref_protein.gff3  -T -o ref_protein.gtf
```

- Count reads using featureCounts.
  - `-Q 0  -M`: Allow multiple-mapping reads that are aligned to different positions. This is good for quantifying duplicated genes.

Count for one sample.
```sh
mkdir -p out_cnt/sep/
featureCounts  -Q 0  -M  -T 8  -s 2  -a ref_protein.gtf  -o out_cnt/sep/S1_Rep1.txt  out_bams/S1_Rep1_fixNH.bam
```

Count for all samples listed in file `list.in_rd`. The result file `out_cnt/joint-cnt.txt` has all samples counted.
```sh
mkdir -p out_cnt/sep/
# cut -f 4 list.in_rd |tail -n +2|perl -e 'while (<>) {chomp;print "featureCounts -Q 0 -M -T 8 -s 2 -a ref_protein.gtf -o out_cnt/sep/$_.cnt out_bams/${_}_fixNH.bam\n";}' > cx2cnt
cut -f 4 list.in_rd|tail -n +2|perl -e 'my @f=<>;chomp(@f);print join(" ", "featureCounts -Q 0 -M -T 8 -s 2 -a ref_protein.gtf -o out_cnt/joint-cnt.txt", map {"out_bams/${_}_fixNH.bam"} @f)."\n";' > cx2cnt
bash cx2cnt
perl -e 'while (<>) {s!out_bams/!!g; s!_fixNH.bam!!g; print;}' out_cnt/joint-cnt.txt > out_cnt/joint-cnt_rename.txt
```
Convert read counts to TPM.
```sh
Rscript cnvt_featureCounts_to_tpm.r out_cnt/joint-cnt_rename.txt out_cnt/joint-tpm.txt
```

## (3) Detect differentially expressed genes with DESeq2.
- Prepare sample meta file `list.sample_meta`.
```sh
echo -e "sample\tgroup" > list.sample_meta; tail -n +2 list.in_rd|deal_table.pl -column 3,0 >> list.sample_meta;
```
Example of `list.sample_meta`.
| sample           | group      |
|------------------|------------|
| S1G1T1P10D\_Rep1 | S1G1T1P10D |
| S1G1T1P10D\_Rep2 | S1G1T1P10D |
| S1G1T1P18D\_Rep1 | S1G1T1P18D |
| S1G1T1P18D\_Rep2 | S1G1T1P18D |


- Prepare a file listing all comparisons `list.comparison`.

|   group1   | group2     | outPrefix                  |
|------------|------------|----------------------------|
| S1G1T1P10D | S1G1T1P18D | PI296341\_flesh\_10D\_18D  |
| S1G1T1P18D | S1G1T1P26D | PI296341\_flesh\_18D\_26D  |


- Compare gene expression for each comparison.
```sh
mkdir -p DEGs/sep/
tail -n +2 list.comparison| perl -e 'while (<>) {chomp; my @a=split; print "Rscript run_deseq2_tpm.r -c joint-cnt_rename.txt  -s list.sample_meta  -t joint-tpm.txt  -g group  --baseline $a[0]  --treatment $a[1]  -o DEGs/sep/res-$a[2]\n";}' > cx3deg
bash cx3deg
```


