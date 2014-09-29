### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

# Basic settings (parameters). 
db='/data/Sunhh/database/db_bwa/Watermelon/WM97_v6.cds.fasta'
db_lab='WM97_cds'

exe_bwa="bwa"
exe_samtools="samtools"
cpuN=10
mem_limit='5G'

para_bwaAln="-t $cpuN -n 4 -o 1 -e 2"
para_bwaSpe=""
para_bwaSse=""
para_bamSort="-@ $cpuN -m $mem_limit"

# Used functions. 
function bwaPE {
	inFq1_tmp=$1
	inFq2_tmp=$2
	oBamPre_tmp=$3
	db_tmp=$4
	db_lab_tmp=$5

	[ -n "$oBamPre_tmp" ] || oBamPre_tmp="bwaOutPre"
	[ -n "$db_tmp" ]      || db_tmp="$db"
	[ -n "$db_lab_tmp" ]  || db_lab_tmp="$db_lab"

	if [[ -f "$inFq1_tmp" && -f "$inFq2_tmp" ]]
	then
		tsmsg "[Rec] Run paired alignment for $inFq1_tmp and $inFq2_tmp to DB=$db_tmp"
	else
		tsmsg "[Err] Some file absent! Check"
		tsmsg "[Err] $inFq1_tmp"
		tsmsg "[Err] $inFq2_tmp"
		exit 1
	fi

	oSai1_tmp="${inFq1_tmp}.${db_lab}.sai"
	oSai2_tmp="${inFq2_tmp}.${db_lab}.sai"

	exe_cmd "$exe_bwa aln $para_bwaAln -f $oSai1_tmp $db_tmp $inFq1_tmp"
	exe_cmd "$exe_bwa aln $para_bwaAln -f $oSai2_tmp $db_tmp $inFq2_tmp"
	exe_cmd "$exe_bwa sampe $para_bwaSpe $db_tmp $oSai1_tmp $oSai2_tmp $inFq1_tmp $inFq2_tmp | $exe_samtools view -bSh -o $oBamPre_tmp.bam -"
	[ -n "$para_bamSort" ] && exe_cmd "$exe_samtools sort $para_bamSort $oBamPre_tmp.bam $oBamPre_tmp.srt"
	[ -n "$para_bamSort" ] && exe_cmd "$exe_samtools index $oBamPre_tmp.srt.bam"

	exe_cmd "rm $oSai1_tmp $oSai2_tmp $oBamPre_tmp.bam"
	unset inFq1_tmp inFq2_tmp oBamPre_tmp db_tmp db_lab_tmp oSai1_tmp oSai2_tmp
}

function bwaSE {
	inFqS_tmp=$1
	oBamPre_tmp=$2
	db_tmp=$3
	db_lab_tmp=$4

	[ -n "$oBamPre_tmp" ] || oBamPre_tmp="bwaOutPre"
	[ -n "$db_tmp" ]      || db_tmp="$db"
	[ -n "$db_lab_tmp" ]  || db_lab_tmp="$db_lab"

	if [[ -f "$inFqS_tmp" ]]
	then
		tsmsg "[Rec] Run single alignment for $inFqS_tmp to DB=$db_tmp"
	else
		tsmsg "[Err] Some file absent! Check"
		tsmsg "[Err] $inFqS_tmp"
		exit 1
	fi

	oSaiS_tmp="${inFqS_tmp}.${db_lab}.sai"

	exe_cmd "$exe_bwa aln $para_bwaAln -f $oSaiS_tmp $db_tmp $inFqS_tmp"
	exe_cmd "$exe_bwa sampe $para_bwaSse $db_tmp $oSaiS_tmp $inFqS_tmp | $exe_samtools view -bSh -o $oBamPre_tmp.bam -"
	[ -n "$para_bamSort" ] && exe_cmd "$exe_samtools sort $para_bamSort $oBamPre_tmp.bam $oBamPre_tmp.srt"
	[ -n "$para_bamSort" ] && exe_cmd "$exe_samtools index $oBamPre_tmp.srt.bam"

	exe_cmd "rm $oSaiS_tmp $oBamPre_tmp.bam"
	unset inFqS_tmp oBamPre_tmp db_tmp db_lab_tmp oSaiS_tmp
}

# Input files and commands. 

bwaPE di_L2_I035.R1.clean.fastq.noRR di_L2_I035.R2.clean.fastq.noRR di_toWM97cds_P
bwaPE gao_L2_I036.R1.clean.fastq.noRR gao_L2_I036.R2.clean.fastq.noRR gao_toWM97cds_P

tsmsg "[Rec] All done."

