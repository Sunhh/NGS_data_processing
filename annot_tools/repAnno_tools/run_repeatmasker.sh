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
pl_buildSum='/data/Sunhh/src/Annot/repeatmasker/RepeatMasker/util/buildSummary.pl'

cpuN=20
repLib="allRepeats_v1.lib"
seqfile="P1All.scf.fa"
minLen=100

tsmsg "[Rec] All start."

# exe_cmd "$exe_rpmk -s -x -lib $repLib $seqfile -nolow -norna -no_is -pa $cpuN -a 1>stdout.RepMsk_unknown 2>stderr.RepMsk_unknown"
tsv_genom='P1Genom_Gt5h.scf.fa.tsv'
exe_cmd "perl $pl_buildSum -useAbsoluteGenomeSize -genome $tsv_genom $seqfile.out > $seqfile.out.NuclGenom.summary"

tsmsg "[Rec] All done."



tsmsg "[Rec] All done."
