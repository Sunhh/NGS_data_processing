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

# Parameter settings for different genomes.
#mkPref="21TJS6_FalT1"
#gnFa="/data/Sunhh/wmhifi/analysis/gene_prediction/in_genomes/21TJS6.fa"
#gnPref="21TJS6"
#fnList="/data/Sunhh/wmhifi/analysis/gene_prediction/07_evm/list.gnPref_mkPref_gnFa_gff"
mkPref=$1
gnFa=$2
#gnPref=$3
fnList=$3
### Check the "Prepare protein alignments" section!!!
### Check if you have "EST alignments" and revise the '-no_est' parameter in PL_pipeEVM!!!

# Parameter settings in general
CONDA_activate_liftoff="source /data/Sunhh/install/anaconda3/bin/activate liftoff"
CONDA_deactivate="conda deactivate"

export EVM_HOME=/data/Sunhh/src/annotation/evm/EVidenceModeler-1.1.1/

PL_trim2cds="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/fmt_gff_trim2CDS.pl"
PL_retGoodGff="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/ret_good_model_from_liftoff_gff3.pl"
PL_fitEVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/fit_evm_pred.pl"
PL_cnvtMk2EVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/cnvt_maker2evmProtGff3.pl"
PL_rmHeadPartCDS="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/rm_head_partial_frame_inGff.pl"
PL_getCDSgff="perl /home/Sunhh/tools/github/NGS_data_processing/temp/get_cds_from_gff3.pl"
# PL_pipeEVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/pipe_revise_byEVM.pl"
PL_pipeEVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/pipe_revise_byEVM.pl -no_est "
PL_get_longCDS_in_gff="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/get_longCDS_in_gff.pl"
PL_runCmd="run_cmd_in_batch.pl"
PL_dealFas="deal_fasta.pl"
PL_dealTab="deal_table.pl"
PL_dealGff="perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl"

EXE_liftoff="liftoff"

### Prepare "param.cfg" for EVM pipeline.
printf "# Configuration file to list nearly all required information for cross-species gene revision process.
### Input section:
in_refList       ref_list
in_addList       add_list
in_parList       param_list

### Software, parameters, and database section:
out_pref         merged
max_olapLen      0
max_olapRat      0

exe_liftoff              $EXE_liftoff
fixPar_liftoff           -copies -polish -exclude_partial -p 40 
pl_get_longCDS_in_gff    $PL_get_longCDS_in_gff
pl_ret_good_model        $PL_retGoodGff
pl_rm_HPartFrame         $PL_rmHeadPartCDS
pl_run_cmd_inBatch       $PL_runCmd
cpuN                     80
##### pl_fix_liftoff_stop      perl /data/Sunhh/cpepo/gene_prediction/x01_cross_species_chk/by_merge_genes/tools/dropped/fix_liftoff_stopCodon.pl
#####  This script is not required since I can use \"-polish\" parameter in the v1.6.3 version.
pl_deal_fasta            $PL_dealFas
pl_deal_table            $PL_dealTab
pl_cds_from_gff3         $PL_getCDSgff
pl_deal_gff3             $PL_dealGff

EVM_HOME                 $EVM_HOME
pl_fitEvmPred            $PL_fitEVM
fn_evm_weight            $PWD/evm_weight.txt
EVM_abPred_weight        6
EVM_transPred_weight     5
fn_evm_protAln           $PWD/prot_aln.evm.gff3
fn_evm_estAln            $PWD/est_aln.evm.gff3
para_evm_partition       --segmentSize 800000  --overlapSize 10000
" > param.cfg

# Prepare genome fasta file
# [[ -e db/ ]] || mkdir db/
# exe_cmd "cp -p ${gnFa} ./db/${gnPref}.fa"
[[ -e tmp/ ]] || mkdir tmp/

# Prepare predicted gene models.
### From maker (include abinit with IPR domains)
# /data/Sunhh/cpepo/analysis/gene_prediction/run_C31/06_maker/maker_results/3fal2/iprscan/comb_maker_noFa.gff3
# GFF_pred="/data/Sunhh/wmhifi/analysis/gene_prediction/06_maker/maker_results/${mkPref}/iprscan/comb_maker_noFa.gff3"

# Run Liftoff to transfer predictions.
# exe_cmd "$CONDA_activate_liftoff"
printf "" > add_list
while IFS=$'\t' read -r -a myArr
do
  if [[ ${myArr[1]} == $mkPref ]];
  then
    exe_cmd "$PL_trim2cds ${myArr[3]} > tmp/${myArr[1]}.trim2cds.gff3"
    echo -e "${myArr[1]}\t${myArr[2]}\t$PWD/tmp/${myArr[1]}.trim2cds.gff3\tNA\tNA" > ref_list
  fi
  exe_cmd "$PL_trim2cds  ${myArr[3]} > tmp/${myArr[1]}.trim2cds.gff3"
  echo -e "${myArr[1]}\t${myArr[2]}\t$PWD/tmp/${myArr[1]}.trim2cds.gff3\tNA\tNA" >> add_list
  # exe_cmd "$EXE_liftoff -polish -copies -sc 0.97 -exclude_partial -s 0.97 -a 0.95 -p 40 -dir tmp/im.${myArr[1]} -o tmp/on.${myArr[1]}.ori.gff3 -u tmp/unmapped.${myArr[1]} -g tmp/${myArr[1]}.trim2cds.gff3 ./db/${gnPref}.fa ${myArr[2]}"
  # exe_cmd "$PL_retGoodGff  tmp/on.${myArr[1]}.ori.gff3_polished > tmp/on.${myArr[1]}.gff3"
done < "$fnList"
# exe_cmd "$CONDA_deactivate"

# Skip evidences from PASA alignments.
### perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/pasa_gff_to_alnGff.pl -notPasa -addTag 'tdcPrj:' ../../db/from_pasa/C31prj_pasa.assemblies.fasta.transdecoder.genome.ID.gff3 | grep -v ^$ > pasaPrj.tdc.gff3
### perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/pasa_gff_to_alnGff.pl -notPasa -addTag 'tdcSra:' ../../db/from_pasa/C31sra_pasa.assemblies.fasta.transdecoder.genome.ID.gff3 | grep -v ^$ > pasaSra.tdc.gff3

# Combine gene predictions.

# Prepare protein alignments.
cat /data/Sunhh/wmhifi/analysis/gene_prediction/03_aln_prot/run_${gnPref}/*.spaln.s2.4maker.gff3 | $PL_cnvtMk2EVM > prot_aln.evm.gff3

# Prepare EST alignments.
### Skip this for now.
# cat ../../db/from_pasa/C31_compreh_prj.4maker.gff3 ../../db/from_pasa/C31_compreh_sra.4maker.gff3 | perl cnvt_maker2evmProtGff3.pl -out_featSource "" -out_featType "EST_match" > est_aln.evm.gff3

# Prepare weight matrix file.
### evm_weight.txt
# echo -e "ABINITIO_PREDICTION\tmaker\t6" > evm_weight.txt
# echo -e "ABINITIO_PREDICTION\tLiftoff\t5" >> evm_weight.txt
#echo -e "OTHER_PREDICTION\ttransdecoder\t3" >> evm_weight.txt
echo -e "PROTEIN\tspliced_protein_alignments\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tblat-DB_pasa\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tgmap-DB_pasa\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tassembler-DB_pasa" >> evm_weight.txt

# Run EVM
exe_cmd "$CONDA_activate_liftoff"
exe_cmd "$PL_pipeEVM -in_cfg ./param.cfg  -bestModel evm -in_refList ref_list -in_addList add_list -out_pref evmMerged"
exe_cmd "$CONDA_deactivate"

