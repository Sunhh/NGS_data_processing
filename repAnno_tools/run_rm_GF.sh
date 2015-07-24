### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

pl_GF="$HOME/tools/github/NGS_data_processing/repAnno_tools/rm_geneFrag.pl"
pl_PE="/data/Sunhh/P1_annot/01.P1_repeat/05.rmGeneFrag/tools/ProtExcluder1.1/ProtExcluder.pl"
pl_dealFa="$HOME//tools/github/NGS_data_processing/deal_fasta.pl"

dbProt='uniprot_sprot_plants_rmTransProt.fa'
dbProt='/data/Sunhh/database/db_fasta/uniprot/20140917/uniprot_sprot_plants_rmTransProt.fa'

inLibLis="inLibLis"
cpuN=30

tsmsg "[Rec] Start"

exe_cmd "perl $pl_GF -evalue 1e-2 -rawLibLis $inLibLis -cpuN $cpuN -dbProt $dbProt -pl_ProtExcluder $pl_PE -pl_dealFa $pl_dealFa"

tsmsg "[Rec] All done."

