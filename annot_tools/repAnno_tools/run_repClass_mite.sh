### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_cmd "/share/app/Annotation/repeatmodeler/RepeatModeler/RepeatClassifier -consensi MITE_rmGen_chop.fa -engine ncbi"
