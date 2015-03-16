### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_rpmk="/data/Sunhh/src/Annot/repeatmasker/RepeatMasker/RepeatMasker"

cpuN=10
repLib="KnownRepeats_v1.lib"
seqfile="P1All.scf.fa"
minLen=100

tsmsg "[Rec] All start."

exe_cmd "$exe_rpmk -x -lib $repLib $seqfile -nolow -norna -no_is -pa $cpuN -a 1>stdout.RepMsk_known 2>stderr.RepMsk_known"

tsmsg "[Rec] All done."
