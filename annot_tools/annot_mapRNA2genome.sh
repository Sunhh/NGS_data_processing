### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

### Settings of EXEs 
dirCuff="/data/Sunhh/src/Annot/Cufflinks/cufflinks-2.2.1.Linux_x86_64"
path_tophat2="tophat2"

mode="single"

### Settings of Personal data
cpuN=60
dbBwt2="db/PG1All_v2_Scf.unmask.fa"
dbTag="toP1G"

oPref="AllToPG1"
para_th="-p $cpuN --library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30"
para_cl="-p $cpuN --min-intron-length 30 --min-frags-per-transfrag 5 "

inDir="reads"
inFileLis="cult_inLis"
if [ "$inDir" != "" ]
then
	inDir="$inDir/"
fi
files=()
while read line
do
	files+=("${inDir}$line")
done <"$inFileLis"


### Running 
# Yi recommanded to merge all reads together and run tophat only one time. 
if [ "$mode" == "pair" ]
then
	exit 1; 
elif [ "$mode" = "single" ]
then
	# Prepare input reads' list 
	inRd1=""
	for tmpID in ${files[@]}
	do
		if [ "$inRd1" != "" ]
		then
			inRd1="${inRd1},$tmpID"
		else
			inRd1="$tmpID"
		fi
	done
	# Running programs. 
	exe_cmd "$path_tophat2 $para_th -o ${oPref}_thout $dbBwt2 $inRd1"
	exe_cmd "$dirCuff/cufflinks $para_cl -o ${oPref}_clout ${oPref}_thout/accepted_hits.bam"
	exe_cmd "cp -p ${oPref}_clout/transcripts.gtf step4.transcripts.gtf"
else 
	exit 2
fi




tsmsg "[Rec] All done."
