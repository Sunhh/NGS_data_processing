inFa=$1
dbFa=$2
cutN=$3
pref=$4

if [ -n "$pref" ] 
then
	deal_fasta.pl $inFa -cut $cutN -cut_prefix $pref
	ls ${pref}_cutted/*.fasta | perl -e ' while (<>) { chomp; print "blastp -query $_ -out $_.blast -db all.fa -evalue 1e-5 -max_target_seqs 1000 -outfmt 6 -num_threads 3 \n"; } ' > cmd_list_blastp_${pref}
	echo "CMD:    nohup run_cmd_in_batch.pl -cpuN 60 cmd_list_blastp_${pref} > scrn.sep_${pref}"
else
	echo "bash $0    in.fa    db_name    cut_num    cut_pref"
	exit
fi

