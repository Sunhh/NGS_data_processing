pl_dropB="/home/Sunhh/tools/github/NGS_data_processing/drop_dup_both_end.pl"

subseqL=0
subseqS=0

# for inPre in NSP306_3hb NSP306_5hb NSP306_1kb NSP306_5kb NSP306_10kb NSP306_15kb
for inPre in NSP306_3hbr2
do
	outPre=$inPre
	inFq1="${inPre}_R1.fastq.gz"
	inFq2="${inPre}_R2.fastq.gz"
	cmd="perl $pl_dropB -opre $outPre $inFq1 $inFq2 -subseq $subseqL -subseqS $subseqS -rcDup"
	echo "[Rec][$(date)] $cmd"
	eval $cmd
done
echo "[Rec][$(date)] Finished." 

# drop_dup_both_end.pl -opre CG_MiSeq02 CG_MiSeq02_R1.fastq.gz CG_MiSeq02_R2.fastq.gz -subseq 100 

