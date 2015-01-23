### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

roundN=4

path_gffMerge=gff3_merge
path_fasMerge=fasta_merge

db_idx=PG1All_v2_Scf.unmsk_master_datastore_index.log

exe_cmd "mkdir ../result/r$roundN"
exe_cmd "$path_gffMerge -d $db_idx -g -o ../result/r${roundN}/r${roundN}_maker.gff3"
exe_cmd "$path_gffMerge -d $db_idx -n -o ../result/r${roundN}/r${roundN}_all.gff3"
exe_cmd "$path_fasMerge -d $db_idx -o ../result/r${roundN}/r${roundN}"

