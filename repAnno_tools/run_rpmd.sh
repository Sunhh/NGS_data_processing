### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}


dir_rpmd='/data/Sunhh/src/Annot/repeatmodeler/RepeatModeler'
exe_rpmd="$dir_rpmd/RepeatModeler"
exe_rpmdBD="$dir_rpmd/BuildDatabase"
exe_rpmk='/data/Sunhh/src/Annot/repeatmasker/RepeatMasker/RepeatMasker'

pl_dealFa='/home/Sunhh/tools/github/NGS_data_processing/deal_fasta.pl'

refFa='P3Genom_Gt5h.scf.fa'
repDb='all_LTR_MITE.lib'

exe_cmd "$exe_rpmk -lib $repDb $refFa -x -nolow -norna -no_is -pa 40 -a 1>stdout.RepMsk 2>stderr.RepMsk"
exe_cmd "perl $pl_dealFa -listSite '[ATGCNatgcn]+' $refFa.masked > $refFa.um_list"
exe_cmd "perl $pl_dealFa $refFa -drawByList -drawList $refFa.um_list -drawLcol 0,2,3 > $refFa.um"

exe_cmd "$exe_rpmdBD -name um -engine ncbi $refFa.um"
exe_cmd "$exe_rpmd -database um -pa 40 "

