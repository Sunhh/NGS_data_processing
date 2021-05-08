### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

pl_goodWrn="$HOME/tools/github/NGS_data_processing/annot_tools/snap_good_wrn_by_valid.pl"
dir_snap='/Data/Sunhh/src/annotation/maker/maker/exe/snap/'
export PATH="$PATH:$dir_snap"

ovl_len=2000
# ln -s ../../../05.run_maker/result/r4/02.fullIntron/r4_finalUse.gff3 ./

in_gff=r5_maker_fullIntron_Island2kb_wiFa.gff3
# woFa_gff=r4_finalUse.gff3
# genom_fa=/Data/Sunhh/seldom/forWushan/rose_genome/gene_prediction/db/samantha_ctg.fasta
# in_gff="${woFa_gff}.wiFa.gff"
org_name=P1Genom
useTag="r5FullIntronIsland2kb"
ohmm_name=${org_name}.${useTag}.hmm

cat $woFa_gff > $in_gff
echo "##FASTA" >> $in_gff
cat $genom_fa >> $in_gff

exe_cmd "maker2zff -n $in_gff"
# exe_cmd "maker2zff -x 0 -o 1 -c 1 -l 30 $in_gff"
exe_cmd "fathom genome.ann genome.dna -validate   > genome.fat_valid"
exe_cmd "fathom genome.ann genome.dna -gene-stats > genome.fat_stats"
exe_cmd "fathom -categorize $ovl_len genome.ann genome.dna"
exe_cmd "perl $pl_goodWrn -valid genome.fat_valid -wrnDna wrn.dna -wrnAnn wrn.ann -yaN 2 -outPref wrnGood"
exe_cmd "cat uni.dna wrnGood.dna > ${useTag}.dna ; cat uni.ann wrnGood.ann > ${useTag}.ann"
exe_cmd "fathom ${useTag}.ann ${useTag}.dna -validate   > ${useTag}.fat_valid"
exe_cmd "fathom ${useTag}.ann ${useTag}.dna -gene-stats > ${useTag}.fat_stats"
exe_cmd "fathom -export $ovl_len -plus ${useTag}.ann ${useTag}.dna"
exe_cmd "fathom export.ann export.dna -validate   > export.fat_valid"
exe_cmd "fathom export.ann export.dna -gene-stats > export.fat_stats"
exe_cmd "mkdir parameters"
exe_cmd "cd parameters ; forge ../export.ann ../export.dna ; cd ../"
exe_cmd "hmm-assembler.pl $org_name parameters > $ohmm_name"

tsmsg "rm -rf genome.??? wrn.??? alt.??? olp.??? err.??? wrnGood.??? parameters/"
