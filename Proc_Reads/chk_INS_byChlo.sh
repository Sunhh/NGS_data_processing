
#db="db/DiSP_slct.scf.fa"
#db_tag="toDiSPslct"
db="db/Cs_chloDNA.fasta"
db_tag="toCsChlo"

pl_samF=/home/Sunhh/tools/github/NGS_data_processing/sam_filter.pl

for inPref in HKC_15_20kb_pTr HKC_8_10kb_pTr HWB_15_20kb_pTr HWB_8_10kb_pTr
do
	rawFq1="${inPref}_R1.fq"
	rawFq2="${inPref}_R2.fq"
	inFq1="$rawFq1.tt"
	inFq2="$rawFq2.tt"
	awk ' NR%40==1 || NR%40==2 || NR%40==3 || NR%40==4 ' $rawFq1 | deal_fastq.pl -frag 0-0 -frag_r -frag_c > $inFq1
	awk ' NR%40==1 || NR%40==2 || NR%40==3 || NR%40==4 ' $rawFq2 | deal_fastq.pl -frag 0-0 -frag_r -frag_c > $inFq2
	# cat $rawFq1 | deal_fastq.pl -frag 0-0 -frag_r -frag_c > $inFq1
	# cat $rawFq2 | deal_fastq.pl -frag 0-0 -frag_r -frag_c > $inFq2
	bwa aln -t 20 -n 4 -o 1 -e 2 $db $inFq1 > $inFq1.$db_tag.sai
	bwa aln -t 20 -n 4 -o 1 -e 2 $db $inFq2 > $inFq2.$db_tag.sai
	bwa sampe -s $db $inFq1.$db_tag.sai $inFq2.$db_tag.sai $inFq1 $inFq2 | samtools view -F 4 -Sh - | perl $pl_samF -bothEnd_pair | samtools view -bSh -o ${inPref}.$db_tag.bam -
	samtools sort ${inPref}.$db_tag.bam ${inPref}.$db_tag.srt
	# samtools index ${inPref}.$db_tag.srt.bam 
	samtools view ${inPref}.$db_tag.srt.bam | sam_filter.pl -h2diff_F | awk ' ($4 < 30000 && $4+$9 > 30000 && $4+$9 < 80000) || ($4 > 30000 && $4+$9 < 30000 && $4 < 80000) { print $9 } ' > ${inPref}.$db_tag.srt.bam.slct.ins
	rm $inFq1.$db_tag.sai $inFq2.$db_tag.sai ${inPref}.$db_tag.bam 
	# samtools view ${inPref}.$db_tag.srt.bam | perl $pl_samF -h2diff_F | cut -f 9 > ${inPref}.$db_tag.srt.bam.ins
done

