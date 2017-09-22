### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

pl_bpTrans=$HOME/tools/github/NGS_data_processing/bp0_2_bp6.pl
pl_slctGff=$HOME/tools/github/NGS_data_processing/annot_tools/slct_maker_gff3.pl

dbPath=prot_db/uniprot_sprot_plants.fas
dbTag="Sprot"
cpuN=50

maxDist2Edge=9

rN=2

inProtFa="r${rN}_maker.prot.fa"
inGff="r${rN}_maker.gff3"
outGff="r${rN}_maker_good.gff3"

exe_cmd "blastp -evalue 1e-10 -query $inProtFa -db $dbPath -num_threads $cpuN -out ${inProtFa}.to${dbTag}.bp0"
exe_cmd "perl $pl_bpTrans -in ${inProtFa}.to${dbTag}.bp0 -out ${inProtFa}.to${dbTag}.bp6"
exe_cmd "awk ' \$7 <= $maxDist2Edge+1 && \$9 <= $maxDist2Edge+1 && \$10 >= \$14-$maxDist2Edge && \$8 >= \$13-$maxDist2Edge && \$3 >= 60 ' ${inProtFa}.to${dbTag}.bp6 > ${inProtFa}.to${dbTag}.bp6.good" 
exe_cmd "perl $pl_slctGff ${inProtFa}.to${dbTag}.bp6.good $inGff > $outGff"


