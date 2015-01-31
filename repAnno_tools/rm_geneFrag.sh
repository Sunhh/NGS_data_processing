### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_cmd "mkdir Final/"
exe_cmd "mkdir MED/"

# for inLib in allLTR.libnPr2 ModelerID.libnPr2 ModelerUnknown.libnPr2 use_MITE.libnPr2
for inLib in allLTR.lib ModelerID.lib ModelerUnknown.lib use_MITE.lib
do
	exe_cmd "blastx -evalue 1e-5 -db uniprot_sprot_plants_rmTransProt.fa -query $inLib -num_threads 40 -out ${inLib}.toPDB.bx0"
	exe_cmd "perl /share/app/Annotation/repeatmodeler/ProtExcluder1.1/ProtExcluder.pl ${inLib}.toPDB.bx0 $inLib"
	exe_cmd "mv ${inLib}noProtFinal Final/"
	exe_cmd "mv ${inLib}.toPDB.bx0* ${inLib}.ssi ${inLib}libnPr MED/"
done

tsmsg "[Rec] All done."; 
