### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dir_aug='/data/Sunhh/src/Annot/maker/2.31.9/maker/exe/augustus/'

export PATH="$dir_aug/bin/:$PATH"


pl_zff2augustus=$HOME/tools/github/NGS_data_processing/annot_tools/zff2augustus_gbk.pl
pl_fixZff2Aug=$HOME/tools/github/NGS_data_processing/annot_tools/fix_1bpLoc_by_zff2Gb.pl
para_rmLis=""
pl_gff2gb=$dir_aug/scripts/gff2gbSmallDNA.pl
pl_randSplit=$dir_aug/scripts/randomSplit.pl
pl_newSpec=$dir_aug/scripts/new_species.pl
pl_copySpec=$HOME/tools/github/NGS_data_processing/annot_tools/copy_species.pl
pl_filtGene=$dir_aug/scripts/filterGenes.pl
pl_optPara=$dir_aug/scripts/optimize_augustus.pl

# pl_rmRed=$HOME/tools/github/NGS_data_processing/annot_tools/keep_nonRedundant_list.pl
pl_rmRedFromProt=$HOME/tools/github/NGS_data_processing/annot_tools/rmRedunt_inputProt.pl
pl_dealTbl=$HOME/tools/github/NGS_data_processing/deal_table.pl
pl_dealFas=$HOME/tools/github/NGS_data_processing/deal_fasta.pl


exe_etrain=etraining
exe_augustus=augustus

# Need export.dna and export.ann files, which come from snap output. 
#   orgRef should be 'generic' for the first run. 
orgRef='BGr2Good'
orgName='BGr3Goods2'
testNum=200
in_genom_fa=export.dna
in_genom_zff=export.ann
in_genom_pep=export.aa
# in_raw_gff=r5_maker_good_noFa.gff3
# in_genom_fa=P1All.scf.fa
in_raw_gb=r2Good.gb
in_gb=train.gb


### Generate list to remove redundancy by proteins similarity
exe_cmd "perl $pl_rmRedFromProt -minIdentity 70 -prot_qry $in_genom_pep -opref ${in_genom_pep}.rmR"
para_rmLis="-rmLis=${in_genom_pep}.rmR.bp6.redund_list"
#exe_cmd "perl $pl_dealFas $in_genom_pep -drawByList -drawList ${in_genom_pep}.rmR.bp6.redund_list -drawLcol 0 -drawWhole -dropMatch > ${in_genom_pep}.rmR.fas"

## Section 1
## Generate a new species in augustus or copy a previous one to the current orgName. Edit information as needed. 
exe_cmd "perl $pl_copySpec --species=$orgName --from_species=$orgRef"
# exe_cmd "perl $pl_copySpec --species=$orgName --from_species=$orgRef --from_trained"
#### Prepare for augustus and etrain : By zff2augustus_gbk.pl Method. 
exe_cmd "perl $pl_zff2augustus $para_rmLis > $in_raw_gb.raw"
#### Prepare for augustus and etrain : By gff2gbSmallDNA.pl Method. 
#exe_cmd "perl $pl_gff2gb $in_raw_gff $in_genom_fa 2000 $in_raw_gb"
#exe_cmd "cp -p ${in_raw_gb} ${in_raw_gb}.raw"

exe_cmd "perl $pl_fixZff2Aug ${in_raw_gb}.raw > $in_raw_gb.fix"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false $in_raw_gb.fix 1>etrain_findErr.std 2>etrain_findErr.err"
exe_cmd "cat etrain_findErr.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb.fix > $in_gb"
exe_cmd "perl $pl_randSplit $in_gb $testNum"
##echo "vim /workdir/laopopo/src/Annot/maker/maker/exe/augustus/config/species/$orgName/${orgName}_parameters.cfg"
##echo "Set stopCodonExcludedFromCDS to false , because our self genbank includes stop codon."
#exit; 

# Section 2
### Training without optimize. If the new species_model comes from a trained species, this step is not requried. 
###   And if the training dataset is bad, say redundant or with some other problem, the further training will make the model worse!!! 
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainRaw.std 2>etrain_trainRaw.err"
# Test the first model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName ${in_gb}.test > goodTrain_dbGood_firsttest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttest.out"
### Till now, I have never successfully trained UTR model, so I give up. 
### Training without optimize. UTR=on
#exe_cmd "$exe_etrain --species=$orgName --UTR=on --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainRawUtr.std 2>etrain_trainRawUtr.err"
# Test the first model by an ab initial prediction
#exe_cmd "$exe_augustus --species=$orgName --UTR=on ${in_gb}.test > goodTrain_dbGood_firsttestUtr.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttestUtr.out"

### Training with optimize. 
# exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
# exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --UTR=on --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
exe_cmd "perl $pl_optPara --cpus=50 --species=$orgName  ${in_gb}.train 1>opt_train.std 2>opt_train.err"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainOpt.std 2>etrain_trainOpt.err"
## Test the second model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName ${in_gb}.test > goodTrain_dbGood_secondtest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_secondtest.out"


#exe_cmd "$exe_augustus --species=$orgName cegma_GoldDb.gff > goodTrain_dbCegma_firsttest.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbCegma_firsttest.out"

