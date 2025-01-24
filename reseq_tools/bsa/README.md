# Pipeline to do BSA analysis

## Call variants into VCF file with GATK.
In this pipeline, no parental sequence data is required.
- Input files: 
 - Read files: `lowBulk_R1.fq.gz`, `lowBulk_R2.fq.gz`, `highBulk_R1.fq.gz`, `highBulk_R2.fq.gz`;
 - Reference genome file: `db/ref_chr.fa` and related BWA indexes.
 - Software: `bwa`, `samtools`, `sentieon` for GATK;
 - Other information: How many individuals were sampled in `low` and `high` bulks for sequencing.
 
- Output files:
 - Filtered VCF file: `two_bulks.vcf.gz`


```sh
mkdir bam/ gvcf/ vcf/;

bwa mem -M -t 30 db/ref_chr.fa  lowBulk_R1.fq.gz  lowBulk_R2.fq.gz  -R "@RG\tID:lowBulk\tSM:lowBulk\tPL:Illumina\tLB:lowBulk\tPU:NA"    | samtools sort -@ 10 -m 10G -o bam/lowBulk.bam
bwa mem -M -t 30 db/ref_chr.fa  highBulk_R1.fq.gz highBulk_R2.fq.gz -R "@RG\tID:highBulk\tSM:highBulk\tPL:Illumina\tLB:highBulk\tPU:NA" | samtools sort -@ 10 -m 10G -o bam/highBulk.bam

samtools index -@ 20 -b bam/lowBulk.bam
samtools index -@ 20 -b bam/highBulk.bam

# steps: 2- Get metrics; 3- Deduplication; 4- BQSR; 5- GVCF;
cat lis1 | perl -e 'use fileSunhh; use LogInforSunhh;
my $exeSen = "/data/shared/sentieon-genomics-202308.03/bin/sentieon";
&runCmd("rm -f cx2_metrics cx3_dedup cx4_bqsr cx5_gvcf");
my $fa = "/data/Sunhh/wmhifi/analysis/snp_related/gatk_snp/db/ref_chr.fa";
while (<>) {chomp;

 my $p="bam/$_";
 my $p1="gvcf/$_";
 &fileSunhh::write2file("cx2_metrics", "$exeSen driver -t 10 -r $fa -i $p.bam --algo AlignmentStat ${p}_stat.txt\n", ">>");
 &fileSunhh::write2file("cx3_dedup", "$exeSen driver -t 100 -i $p.bam --algo LocusCollector --fun score_info $p.SCORE\n", ">>");
 &fileSunhh::write2file("cx3_dedup", "$exeSen driver -t 100 -i $p.bam --algo Dedup --rmdup --score_info $p.SCORE --metrics $p.dedup-metric.txt $p.dedup.bam\n", ">>");
 &fileSunhh::write2file("cx4_bqsr", "$exeSen driver -t 100 -r $fa -i $p.dedup.bam --algo QualCal $p.recal.table\n", ">>");
 &fileSunhh::write2file("cx5_gvcf", "$exeSen driver -t 30 -r $fa -i $p.dedup.bam --algo Haplotyper --emit_mode gvcf $p1.gvcf.gz\n", ">>");
}'


# step 6- Combine GVCF;
cat lis2 | perl -e 'use fileSunhh; use LogInforSunhh;
my $exeSen = "/data/shared/sentieon-genomics-202308.03/bin/sentieon";
&runCmd("rm -f cx6_call");
my $fa = "/data/Sunhh/wmhifi/analysis/snp_related/gatk_snp/db/ref_chr.fa";
my @x1 = <>; chomp(@x1);
my @fG = map {"gvcf/$_.gvcf.gz"} @x1;
&fileSunhh::write2file("cx6_call", "$exeSen driver -t 90 -r $fa --algo GVCFtyper vcf/noRecal.vcf.gz @fG\n", ">>");
&fileSunhh::write2file("cx6_call", "$exeSen driver -t 90 -r $fa --algo GVCFtyper --emit_conf 1 vcf/noRecalEmQ1.vcf.gz @fG\n", ">>");
'

# step 7- hard filter;
mv vcf/noRecal.vcf.gz vcf/emQ30.vcf.gz
bcftools index --tbi --threads 90 vcf/emQ30.vcf.gz

docker run --rm -v $PWD:/mnt broadinstitute/gatk:4.6.1.0 gatk SelectVariants -R /mnt/db/ref_chr.fa -V /mnt/vcf/emQ30.vcf.gz --select-type-to-include SNP -O /mnt/vcf/emQ30-gatkSNP.vcf.gz
docker run --rm -v $PWD:/mnt broadinstitute/gatk:4.6.1.0 rm -f /mnt/vcf/emQ30-gatkSNP.vcf.gz.tbi
bcftools index --tbi --threads 90 vcf/emQ30-gatkSNP.vcf.gz
bcftools view -e "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AC == AN || AC == 0" vcf/emQ30-gatkSNP.vcf.gz --threads 90 | bcftools +fill-tags --threads 90 -Oz -o vcf/emQ30-gatkSNPhard.vcf.gz
bcftools index --tbi --threads 90 vcf/emQ30-gatkSNPhard.vcf.gz
docker run --rm -v $PWD:/mnt broadinstitute/gatk:4.6.1.0 gatk SelectVariants -R /mnt/db/ref_chr.fa -V /mnt/vcf/emQ30.vcf.gz --select-type-to-include INDEL -O /mnt/vcf/emQ30-gatkInDel.vcf.gz
docker run --rm -v $PWD:/mnt broadinstitute/gatk:4.6.1.0 rm -f /mnt/vcf/emQ30-gatkInDel.vcf.gz.tbi
bcftools index --tbi --threads 90 vcf/emQ30-gatkInDel.vcf.gz
bcftools view -e "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || AC == AN || AC == 0" vcf/emQ30-gatkInDel.vcf.gz --threads 90 | bcftools +fill-tags --threads 90 -Oz -o vcf/emQ30-gatkInDelhard.vcf.gz
bcftools index --tbi --threads 90 vcf/emQ30-gatkInDelhard.vcf.gz

bcftools concat -a -Oz -o vcf/emQ30-gatkAllhard.vcf.gz --threads 90 vcf/emQ30-gatkSNPhard.vcf.gz vcf/emQ30-gatkInDelhard.vcf.gz
bcftools index --tbi --threads 90 vcf/emQ30-gatkAllhard.vcf.gz

# Futher filtering:
### Used rules: bi-allelic, not on CLV01_Chr00/Ct/Mt
### Not used rules: No neighboring variants;
bcftools view -t ^CLV01_Chr00,CLV01_ChrMt,CLV01_ChrCt -c2 -C2 vcf/emQ30-gatkAllhard.vcf.gz -Oz -o vcf/emQ30-gatkAllhard-use.vcf.gz
bcftools index --tbi --threads 90 vcf/emQ30-gatkAllhard-use.vcf.gz

```


## Perform BSA with R/QTLseqr;
- Prepare table file from VCF for R/QTLseqr.
 - Keep only high and low bulks.
 - Output file: `vcf/use-2bulk-chrV.table`
```sh
bcftools view vcf/emQ30-gatkAllhard-use.vcf.gz -s lowBulk,highBulk -Ov -o vcf/use-2bulk.vcf
docker run --rm -v $PWD:/mnt broadinstitute/gatk:4.6.1.0 gatk VariantsToTable -R /mnt/db/ref_chr.fa -V /mnt/vcf/use-2bulk.vcf \
 -F CHROM -F POS -F REF -F ALT \
 -GF AD -GF DP -GF GQ -GF PL \
 -O /mnt/vcf/use-2bulk.table
perl -e 'while (<>) { s!^CLV01_Chr0*!!; s!^\t!0\t!; print;}' vcf/use-2bulk.table > vcf/use-2bulk-chrV.table
```

- Perform BSA.
 - Output files: `output-ana-*`; Tables can be put into `example_data/template-QTLseqr_result.xlsx`;
```sh
Rscript run_QTLseqr.r -i vcf/use-2bulk-chrV.table --high_bulk highBulk --low_bulk lowBulk --indvN_high 20 --indvN_low 20 --window_size 2000000 --window_name 2M --plot_chr 6 -o output
```

