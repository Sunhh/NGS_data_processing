### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}


exe_reapr='/home/Sunhh/src/Assemble/REAPR/Reapr_1.0.17/reapr'
cpuN=10
rd_ident=0.99

in_fa='NSP306_Pla03s01GC_Gt5h.scf.fa'
use_faPref='NSP306_Pla03s01GC_Gt5h'
odir="${use_faPref}_15kb"

#long_fq1=/home/Sunhh/Assembly/Sweetpotato/NSP306/FinalFq/NSP306_15kb_rc_sub100_R1.fq
#long_fq2=/home/Sunhh/Assembly/Sweetpotato/NSP306/FinalFq/NSP306_15kb_rc_sub100_R2.fq
long_fq1=/home/Sunhh/Assembly/Sweetpotato/NSP306/FinalFq/NSP306_15kb_pTr_rc_R1.fq
long_fq2=/home/Sunhh/Assembly/Sweetpotato/NSP306/FinalFq/NSP306_15kb_pTr_rc_R2.fq
fq_label=$odir

shrtBamPref=""
longBam="${fq_label}_long.bam"

para_smalt="-n $cpuN -y $rd_ident"
para_pipe="-break a=1 -score f=15"

tsmsg "[Rec] Begin."

cmd="$exe_reapr facheck $in_fa $use_faPref"
exe_cmd $cmd

cmd="$exe_reapr smaltmap $para_smalt ${use_faPref}.fa $long_fq1 $long_fq2 $longBam"
exe_cmd $cmd

cmd="$exe_reapr pipeline $para_pipe ${use_faPref}.fa $longBam $odir $shrtBamPref"
exe_cmd $cmd

tsmsg "[Rec] All done."
