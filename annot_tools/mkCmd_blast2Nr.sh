### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dbFa='/share/nas2/xigua/sunhonghe/database/db_blast/ncbi/nr'
dbTag='toNr'

qryFa=$1
cpuN=10

echo "blastp -outfmt 11 -db $dbFa -evalue 1e-3 -num_alignments 20 -seg yes -num_threads $cpuN -query $qryFa -out $qryFa.$dbTag.asn.1"

# blastp -outfmt 5 -db nr -evalue 1e-3 -num_alignments 20 -seg yes -num_threads 20 -query wcgP_cutted/wcgP_00033.fasta -out wcgP_cutted/wcgP_00033.fasta.xml


