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
outDir="by_minimap2_nucdiff/output/"
refFaFn=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/final_C39.chr.fa
refGffFn="/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/final_C39.chr.gff3"
qryDir=/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/
PL_affGene="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/get_sv_affected_genes.pl"
PL_rmGap="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/remove_gap_var.pl"
PY_sam2delta="python /data/Sunhh/cpepo/analysis/SV_detection/tools/sam2delta.py"
PL_gff2bed="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/gff2bed.pl"
PL_slctVAR="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/select_var.pl"
PL_getComm="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/get_qry_ref_shared_var_nucdiff.pl"
PL_addAHRD="perl /data/Sunhh/cpepo/analysis/SV_detection/tools/add_geneAHRD.pl"
DIR_db="/data/Sunhh/cpepo/analysis/SV_detection/in_genomes/"

[[ -e $outDir ]] || mkdir -p $outDir
cd $outDir/
for cid in C31 C38
do
  # exe_cmd "minimap2 -ax asm20 -t 20 --cs $refFaFn  $qryDir/final_${cid}.chr.fa  > ${cid}.asm20.sam"
  # exe_cmd "$PY_sam2delta ${cid}.asm20.sam  --ref_fa $refFaFn  --qry_fa $qryDir/final_${cid}.chr.fa"
  # exe_cmd "nucdiff --proc 10 --delta_file ${cid}.asm20.sam.delta  $refFaFn $qryDir/final_${cid}.chr.fa ./ ${cid}.ndf"
  cd results/
  # exe_cmd "$PL_gff2bed ${cid}.ndf_ref_snps.gff > ${cid}.ndf_ref_snps.bed"
  # exe_cmd "$PL_gff2bed ${cid}.ndf_ref_struct.gff > ${cid}.ndf_ref_struct.bed"
  # exe_cmd "$PL_gff2bed ${cid}.ndf_query_snps.gff > ${cid}.ndf_query_snps.bed"
  # exe_cmd "$PL_gff2bed ${cid}.ndf_query_struct.gff > ${cid}.ndf_query_struct.bed"
  # exe_cmd "$PL_rmGap 100 ${cid}.ndf_ref_snps.gff      > ${cid}.ndf_ref_snps.rmGap.gff"
  # exe_cmd "$PL_rmGap 100 ${cid}.ndf_ref_struct.gff    > ${cid}.ndf_ref_struct.rmGap.gff"
  # exe_cmd "$PL_rmGap 100 ${cid}.ndf_query_snps.gff    > ${cid}.ndf_query_snps.rmGap.gff"
  # exe_cmd "$PL_rmGap 100 ${cid}.ndf_query_struct.gff > ${cid}.ndf_query_struct.rmGap.gff"

  # exe_cmd "$PL_getComm -qry_tbl ${cid}.ndf_query_snps.rmGap.gff    -ref_tbl ${cid}.ndf_ref_snps.rmGap.gff    -out_prefix ${cid}.ndf_snps_rmGap"
  # exe_cmd "$PL_getComm -qry_tbl ${cid}.ndf_query_struct.rmGap.gff  -ref_tbl ${cid}.ndf_ref_struct.rmGap.gff  -out_prefix ${cid}.ndf_stru_rmGap"
  # exe_cmd "$PL_affGene -refGff  $refGffFn                      -qryGff ${cid}.ndf_snps_rmGap.ref  -upDist 2000 > ${cid}.ndf_snps_rmGap.ref.effG"
  # exe_cmd "$PL_affGene -refGff  $refGffFn                      -qryGff ${cid}.ndf_stru_rmGap.ref  -upDist 2000 > ${cid}.ndf_stru_rmGap.ref.effG"
  # exe_cmd "$PL_affGene -refGff  $qryDir/final_${cid}.chr.gff3  -qryGff ${cid}.ndf_snps_rmGap.qry  -upDist 2000 > ${cid}.ndf_snps_rmGap.qry.effG"
  # exe_cmd "$PL_affGene -refGff  $qryDir/final_${cid}.chr.gff3  -qryGff ${cid}.ndf_stru_rmGap.qry  -upDist 2000 > ${cid}.ndf_stru_rmGap.qry.effG"
  # exe_cmd "$PL_slctVAR  20  ${cid}.ndf_snps_rmGap.ref.effG  > ${cid}.ndf_snps_rmGap.ref.effG.slct1"
  # exe_cmd "$PL_slctVAR  20  ${cid}.ndf_stru_rmGap.ref.effG  > ${cid}.ndf_stru_rmGap.ref.effG.slct1"
  # exe_cmd "$PL_slctVAR  20  ${cid}.ndf_snps_rmGap.qry.effG  > ${cid}.ndf_snps_rmGap.qry.effG.slct1"
  # exe_cmd "$PL_slctVAR  20  ${cid}.ndf_stru_rmGap.qry.effG  > ${cid}.ndf_stru_rmGap.qry.effG.slct1"

  exe_cmd "$PL_addAHRD  $DIR_db/final_C39.ahrd     ${cid}.ndf_snps_rmGap.ref.effG.slct1  >  ${cid}.ndf_snps_rmGap.ref.effG.slct1_ahrd"
  exe_cmd "$PL_addAHRD  $DIR_db/final_C39.ahrd     ${cid}.ndf_stru_rmGap.ref.effG.slct1  >  ${cid}.ndf_stru_rmGap.ref.effG.slct1_ahrd"
  exe_cmd "$PL_addAHRD  $DIR_db/final_${cid}.ahrd  ${cid}.ndf_snps_rmGap.qry.effG.slct1  >  ${cid}.ndf_snps_rmGap.qry.effG.slct1_ahrd"
  exe_cmd "$PL_addAHRD  $DIR_db/final_${cid}.ahrd  ${cid}.ndf_stru_rmGap.qry.effG.slct1  >  ${cid}.ndf_stru_rmGap.qry.effG.slct1_ahrd"

  exe_cmd "deal_table.pl ${cid}.ndf_snps_rmGap.ref.effG.slct1_ahrd -column 0,1-7 | ColLink.pl -f1 ${cid}.ndf_snps_rmGap.qry.effG.slct1_ahrd -keyC1 0 -keyC2 0 -add -Col1 1-7 | deal_table.pl -column 0,6,7,1-3,8-10,4,5,11,12 > ${cid}.ndf_snps_rmGap.slct1.comb"
  exe_cmd "deal_table.pl ${cid}.ndf_stru_rmGap.ref.effG.slct1_ahrd -column 0,1-7 | ColLink.pl -f1 ${cid}.ndf_stru_rmGap.qry.effG.slct1_ahrd -keyC1 0 -keyC2 0 -add -Col1 1-7 | deal_table.pl -column 0,6,7,1-3,8-10,4,5,11,12 > ${cid}.ndf_stru_rmGap.slct1.comb"


  cd -
done
cd -


