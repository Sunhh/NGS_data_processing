function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_java="java"
exe_jar="/data/Sunhh/src/Assemble/Trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar"

cpuN=30
adp_fas="/data/Sunhh/src/Assemble/Trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa"
minLen=25

para_jar="-threads $cpuN"
para_PE="ILLUMINACLIP:$adp_fas:2:30:10:1 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:$minLen"

tsmsg "[Rec] All start."
for inPref in HKC_15_20kb HKC_8_10kb HWB_15_20kb HWB_8_10kb
do
	tsmsg "[Rec] Dealing with $inPref"
	inFq1="${inPref}.p1"
	inFq2="${inPref}.p2"
	oPref="${inPref}"
	logFile="log.${oPref}"
	# para_jarAdd="-trimlog $logFile"
	para_jarAdd=""
	cmd="$exe_java -jar $exe_jar PE $para_jar $para_jarAdd $inFq1 $inFq2 ${oPref}_pTr_R1.fq ${oPref}_sTr_R1.fq ${oPref}_pTr_R2.fq ${oPref}_sTr_R2.fq $para_PE"
	exe_cmd "$cmd"
done

tsmsg "[Rec] All done."


