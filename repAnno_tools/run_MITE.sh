### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

pl_dealFa=$HOME/tools/github/NGS_data_processing/deal_fasta.pl
pl_addClass=$HOME/tools/github/NGS_data_processing/repAnno_tools/add_repClass.pl

pl_MITE=$HOME/src/Annotation/MITE_Hunter/MITE_Hunter_blast216/MITE_Hunter_manager.pl

refFa='P1All.scf.fa'
outG='P1AllGt5hScf'
cpuN=10
grpN=10

tsmsg "Start."

# exe_cmd "mkdir running/ Step8/"

[ -d "running" ] || exe_cmd "mkdir running"
[ -d "Step8" ] || exe_cmd "mkdir Step8"

cd running
ln -s ../$refFa .

exe_cmd "perl $pl_MITE -c $cpuN -n $grpN -S 12345678 -i $refFa -g $outG"
exe_cmd "cp -p *_Step8*.fa ../Step8/"
exe_cmd "cat ${outG}_Step8*.fa > ../MITE_raw.lib"
exe_cmd "perl $pl_dealFa ../MITE_raw.lib -frag_head -frag_width 80 -frag 0-0 | perl $pl_dealFa -chopKey ':\\d+-\\d+\$' | perl $pl_addClass MITE > ../MITE_named.lib"

cd ../

tsmsg "All done."
