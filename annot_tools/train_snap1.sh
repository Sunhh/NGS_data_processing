### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

# Example command: bash train_snap1.sh in_genome.fa  in_genome.pep.gff3  taxID roundID

# Specify for each server. 
dir_NGS=$HOME/tools/github/NGS_data_processing/
dir_snap='/data/Sunhh/src/annotation/maker/maker/exe/snap'
exe_maker2zff='/data/Sunhh/src/annotation/maker/maker/bin/maker2zff'

# Specify for each run
ovl_len=2000
genom_fa=$1 # in_genome.fa
woFa_gff=$2 # r4_finalUse.gff3
org_name=$3 # P1Genom
useTag=$4   # r5FullIntronIsland2kb

### No need to change. 
ohmm_name=${org_name}.${useTag}.hmm

in_gff="wiFa.gff" # r5_maker_fullIntron_Island2kb_wiFa.gff3
cat $woFa_gff > $in_gff
echo "##FASTA" >> $in_gff
cat $genom_fa >> $in_gff

# Paths 
export PATH="$PATH:$dir_snap"
pl_goodWrn="$dir_NGS/annot_tools/snap_good_wrn_by_valid.pl"

exe_cmd "$exe_maker2zff -n $in_gff"
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
