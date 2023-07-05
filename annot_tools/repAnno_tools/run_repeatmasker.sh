### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_rpmk="/data/Sunhh/src/annotation/repeatmasker/RepeatMasker/RepeatMasker"
pl_buildSum="/data/Sunhh/src/annotation/repeatmasker/RepeatMasker/util/buildSummary.pl"

cpuN=40
repLib="D202306.TElib.fa"
seqfile="hap1.fa"
tsv_genom="$seqfile.tsv"

tsmsg "[Rec] All start."

deal_fasta.pl -baseCount $seqfile | awk 'NR > 1 {print $1"\t"$8-$6}' > $tsv_genom

exe_cmd "$exe_rpmk -s -x -lib $repLib $seqfile -nolow -norna -no_is -pa $cpuN -a 1>$seqfile.stdout_RepMsk 2>$seqfile.stderr_RepMsk"
exe_cmd "perl $pl_buildSum -useAbsoluteGenomeSize -genome $tsv_genom $seqfile.out > $seqfile.out.NuclGenom.summary"

tsmsg "[Rec] All done."

