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

outDir="by_minimap2_paftools/output/"
refFaFn=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/final_C39.chr.fa
qryDir=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/
[[ -e $outDir ]] || mkdir -p $outDir
cd $outDir/
for cid in C31 C38
do
  exe_cmd "minimap2 -cx asm20 -t 20 --cs $refFaFn  $qryDir/final_${cid}.chr.fa  > ${cid}.asm20.paf"
  exe_cmd "sort -k6,6 -k8,8n ${cid}.asm20.paf | paftools.js call -L 10000 -s $cid -f $refFaFn - > ${cid}.asm20.vcf"
done
cd -

