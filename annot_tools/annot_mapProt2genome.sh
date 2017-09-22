### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

### Basic settings
# First you need to add splan_path to $PATH
#                   set $ALN_TAB to /path/to/spaln_dir/table
#                   set $ALN_DBS to /path/to/spaln_dir/seqdb
# Second should use genome.mfa file to build database
#  cd $ALN_DBS
#  ./makeidx.pl -ip genome genome.mfa
# Here 'genome' is the name to be used in future alignment. 
# Third (optional), we need 'gff3_merge' to convert all pep_aln.gff3 to the final one. 

args=("${@}")
myBasename=${0##*/}
path2Me=${0}


### Setting paths. 
cpuN=50
pl_cnvt2maker=$HOME/tools/github/NGS_data_processing/annot_tools/cnvt_spaln2makerAln_prot_gff3.pl
path_gffMerge='gff3_merge'

### Personal files. 
inFileLis="qryProtFasLis"
dbID=PG1All_v2
outFmtGff=all_spaln2genom_prot_wiM4.fmt.gff3

[ -n "${args[0]}" ] && inFileLis="${args[0]}"
[ -n "${args[1]}" ] && dbID="${args[1]}"

### Parsing $inFileLis
if [ "$inFileLis" = "" ] || [ ! -f "$inFileLis" ]
then
	echo "########################################################"
	echo "Is the inFileLis [$inFileLis] OK?"
	echo "Try $path2Me qryProtFasLis spaln_dbID"
	echo ""
	echo "########################################################"
	exit 1; 
fi
qryFiles=( $( cat $inFileLis ) )

for qryFa in ${qryFiles[@]}
do
	qryPara="$qryFa"
	exe_cmd "spaln -t$cpuN -M4 -Q7 -O0 -LS -ya012 -o ${qryFa}.spaln.gff3 -d${dbID} $qryPara"
	exe_cmd "perl $pl_cnvt2maker ${qryFa}.spaln.gff3 > ${qryFa}.med.gff3"
	medGff+=("${qryFa}.med.gff3")
done

medGffLine="${medGff[@]}"
exe_cmd "gff3_merge -l -o $outFmtGff $medGffLine"

