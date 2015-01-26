### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dir_snap='/workdir/laopopo/src/Annot/maker/maker/exe/snap'
export PATH="$PATH:$dir_snap"
# in_gff=r1_maker.gff3
in_gff=r1_maker_good.gff3
org_name=SPG
ohmm_name=SPG.snap1Good.hmm

exe_cmd "maker2zff -x 0 $in_gff"
exe_cmd "fathom genome.ann genome.dna -validate   > genome.fat_valid"
exe_cmd "fathom genome.ann genome.dna -gene-stats > genome.fat_stats"
exe_cmd "fathom -categorize 1000 genome.ann genome.dna"
exe_cmd "fathom uni.ann uni.dna -validate   > uni.fat_valid"
exe_cmd "fathom uni.ann uni.dna -gene-stats > uni.fat_stats"
exe_cmd "fathom -export 1000 -plus uni.ann uni.dna"
exe_cmd "fathom export.ann export.dna -validate   > export.fat_valid"
exe_cmd "fathom export.ann export.dna -gene-stats > export.fat_stats"
exe_cmd "mkdir parameters"
exe_cmd "cd parameters ; forge ../export.ann ../export.dna ; cd ../"
exe_cmd "hmm-assembler.pl $org_name parameters > $ohmm_name"


