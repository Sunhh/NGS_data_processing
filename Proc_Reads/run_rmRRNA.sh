#!/bin/sh

### Change this file name each time. 
fn_pref='pref_list'
printCmd=1


function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][CMD_done] $1"
}

function tsmsg {
	echo "[$(date)]$1"
}

function load_list {
	prefs=() 
	while IFS=$'\t' read -r -a myArray
	do
		prefs+=("${myArray[0]}")
	done < "$1"
}

exe_Bwt="bowtie"
para_Bwt="-p 20 -k 1 -v 3"

db_rrna="/share1/db_bowtie/rRNA_silva111"

exe_pl='perl'
#pl_runBwt="~/tools/github/NGS_data_processing/Proc_Reads/run_bowtie.pl"
pl_extractFq="~/tools/github/NGS_data_processing/extract_fq_by_list.pl"

### Read in a file : Start 
prefs=()
if [[ -n $fn_pref ]]
then
	while IFS=$'\t' read -r -a myArray
	do
		prefs+=("${myArray[0]}")
	done < "$fn_pref"
fi
### Read in a file : End 

### Run bowtie : 
for pref in ${prefs[@]}
do
	inFq1="${pref}_sTr_R1.fq"
	oSam1="${pref}_sTr_R1.rrna.aln"
	outFq1="${pref}_rmRRNA.fq"
	if [[ "$printCmd" -eq "1" ]]
	then
		echo "$exe_Bwt   $para_Bwt   $db_rrna   $inFq1   $oSam1"
		echo "$exe_pl   $pl_extractFq   -refLis $oSam1   -srcFq $inFq1   -outFq $outFq1   -mode drop -rdkey   2>> err.extract"
	else
		exe_cmd "$exe_Bwt   $para_Bwt   $db_rrna   $inFq1   $oSam1"
		exe_cmd "$exe_pl   $pl_extractFq   -refLis $oSam1   -srcFq $inFq1   -outFq $outFq1   -mode drop -rdkey   2>> err.extract"
		tsmsg "[Msg] Finish ${pref}"
	fi
done


