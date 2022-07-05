### Basic functions.
function exe_cmd {
  echo "[$(date)][CMD] $1"
  eval "$1"
  if [[ $? -eq "0" ]]
  then
    echo "[$(date)][CMD_done] $1"
  else
    echo "[$(date)][CMD_err] $1"
  fi
}

function tsmsg {
  echo "[$(date)][Msg] $1"
}


export PATH="/data/Sunhh/src/Align/mummer/install/mummer4/bin/:$PATH"
refFaFn=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/final_C39.chr.fa
qryDir=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/
outDir="by_nucdiff/output/"
[[ -e $outDir ]] || mkdir -p $outDir
cd $outDir/
for cid in C31 C38
do
  exe_cmd "nucdiff --proc 10 --nucmer_opt ' -t 20 --batch 1 ' $refFaFn $qryDir/final_${cid}.chr.fa ./ ${cid}.ndf"
done
cd -

