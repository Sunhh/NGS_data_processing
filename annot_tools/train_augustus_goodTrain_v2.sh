### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

# Specify for each server
dir_NGS=$HOME/tools/github/NGS_data_processing/
dir_aug='/data/Sunhh/src/annotation/augustus/Augustus-3.4.0/'
export AUGUSTUS_CONFIG_PATH=${dir_aug}/config/
dir_braker=/Data/Sunhh/src/annotation/braker/BRAKER_v2.0.4/
exe_maker2zff=/data/Sunhh/src/annotation/maker/maker/bin/maker2zff
exe_fathom=/data/Sunhh/src/annotation/maker/maker/exe/SNAP/fathom

# Specify for each run
doRmRed='y'   # If we want to remove redundancy in proteins. 
doOpt='y' # If we want to optimize the parameters.
in_type=$1    # "export|gff". export for snap output, gff for gff3 file input. 
orgRef=$2  # generic or some other old species model ID. 
orgName=$3 # New model ID. 
testNum=$4 # 2000 / 1000 
ovl_len=2000

if [ $in_type == "gff" ]; then
  in_raw_gff=$5  # r5_maker_good_noFa.gff3
  in_genom_fa=$6 # P1All.scf.fa
elif [ $in_type == "export" ]; then
  # Need export.dna and export.ann files, which come from snap output. 
  A=A; 
elif [ $in_type == "gb" ]; then
  echo "I don't accept gb in_type yet."
  exit; 
else
  echo "Failed to parse in_type value [$in_type]" 
  exit; 
fi


# No need to change
dir_augBin=${dir_aug}/bin/
export PATH="$dir_augBin:$PATH"

# Paths
pl_zff2augustus=$dir_NGS/annot_tools/zff2augustus_gbk.pl
pl_fixZff2Aug=$dir_NGS/annot_tools/fix_1bpLoc_by_zff2Gb.pl
pl_setStopCod=$dir_NGS/annot_tools/set_stopCodonFreq.pl
pl_accuracy=$dir_NGS/annot_tools/augustus.accuracy_calculator.pl
pl_gff2cds=$dir_NGS/temp/get_cds_from_gff3.pl

pl_gff2gb=$dir_aug/scripts/gff2gbSmallDNA.pl
pl_randSplit=$dir_aug/scripts/randomSplit.pl
pl_newSpec=$dir_aug/scripts/new_species.pl
pl_copySpec=$dir_NGS/annot_tools/copy_species.pl
pl_filtGene=$dir_aug/scripts/filterGenes.pl
pl_optPara=$dir_aug/scripts/optimize_augustus.pl

# pl_rmRed=$dir_NGS/annot_tools/keep_nonRedundant_list.pl
pl_rmRedFromProt=$dir_NGS/annot_tools/rmRedunt_inputProt.pl
pl_dealTbl=$dir_NGS/deal_table.pl
pl_dealFas=$dir_NGS/deal_fasta.pl

exe_etrain=$dir_augBin/etraining
exe_augustus=$dir_augBin/augustus

# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/gbSmallDNA2gff.pl
# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/gff2gbSmallDNA.pl
# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/getAnnoFasta.pl
# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/gffGetmRNA.pl
# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/gtf2gff.pl
# /data/Sunhh/src/Annotation/maker/maker.2.31/exe/augustus/scripts/pslMap.pl
# /home/Sunhh/tools/github/NGS_data_processing/temp/get_cds_from_gff3.pl

fn_orgParaCfg=$AUGUSTUS_CONFIG_PATH/species/$orgName/${orgName}_parameters.cfg
in_raw_gb="genes_for_${orgName}.gb"
in_gb=train.gb

# Prepare data for training. 
### I'd like to change all data into export data format. 
if [ $in_type == "gff" ]; then
  in_wiFa_gff="${in_raw_gff}.wiFa.gff"
  cat $in_raw_gff > $in_wiFa_gff
  echo "##FASTA" >> $in_wiFa_gff
  cat $in_genom_fa >> $in_wiFa_gff
  exe_cmd "$exe_maker2zff -n $in_wiFa_gff"
  exe_cmd "$exe_fathom -export $ovl_len -plus genome.ann genome.dna"
fi
### Remove redundancy
para_rmLis=""
##### Generate list to remove redundancy by proteins similarity
if [ $doRmRed == "y" ]; then
  exe_cmd "perl $pl_rmRedFromProt -minIdentity 70 -prot_qry export.aa -opref export.aa.rmR"
  para_rmLis="-rmLis=export.aa.rmR.bp6.redund_list"
  # exe_cmd "perl $pl_dealFas export.aa -drawByList -drawList export.aa.rmR.bp6.redund_list -drawLcol 0 -drawWhole -dropMatch > export.aa.rmR.fas"
fi
exe_cmd "perl $pl_zff2augustus $para_rmLis > ${in_raw_gb}.raw"

# Section 1. Prepare data. 
### Prepare species configuration in augustus or copy a previous one to the current orgName. 
if [ $orgRef == "generic" ]; then
  exe_cmd "perl $pl_copySpec --species=$orgName --from_species=generic"
else
  exe_cmd "perl $pl_copySpec --species=$orgName --from_species=$orgRef --from_trained"
fi
#### Prepare for augustus and etrain : By zff2augustus_gbk.pl Method. 
exe_cmd "perl $pl_fixZff2Aug ${in_raw_gb}.raw > $in_raw_gb.fix"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} $in_raw_gb.fix 1>etrain_findErr.std 2>etrain_findErr.err"
exe_cmd "cat etrain_findErr.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb.fix > $in_gb"
exe_cmd "perl $pl_randSplit $in_gb $testNum"
exe_cmd "perl $pl_setStopCod etrain_findErr.std $fn_orgParaCfg"
##echo "Set stopCodonExcludedFromCDS to false , because our self genbank includes stop codon."

# Section 2. Train for the first time. 
### Training without optimize. If the new species_model comes from a trained species, this step is not requried. 
###   And if the training dataset is bad, say redundant or with some other problem, the further training will make the model worse!!! 
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} ${in_gb}.train 1>etrain_trainRaw.std 2>etrain_trainRaw.err"
# Test the first model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName  --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH}  ${in_gb}.test > train_firsttest.out"
exe_cmd "grep -A 22 Evaluation train_firsttest.out"
exe_cmd "perl $pl_accuracy train_firsttest.out"
### Till now, I have never successfully trained UTR model, so I give up. 
### Training without optimize. UTR=on
#exe_cmd "$exe_etrain --species=$orgName --UTR=on --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainRawUtr.std 2>etrain_trainRawUtr.err"
# Test the first model by an ab initial prediction
#exe_cmd "$exe_augustus --species=$orgName --UTR=on ${in_gb}.test > goodTrain_dbGood_firsttestUtr.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttestUtr.out"

# Section 3. Optimization. 
if [ $doOpt == "y" ]; then
  exe_cmd "perl $pl_optPara --cpus=100 --rounds=5 --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} --species=$orgName  ${in_gb}.train  1>opt_train.std 2>opt_train.err"
  # exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
  # exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --UTR=on --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
  exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} ${in_gb}.train 1>etrain_trainOpt.std 2>etrain_trainOpt.err"
  ## Test the second model by an ab initial prediction
  exe_cmd "$exe_augustus --species=$orgName  --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH}  ${in_gb}.test > train_secondtest.out"
  exe_cmd "grep -A 22 Evaluation train_secondtest.out"
  exe_cmd "perl $pl_accuracy train_secondtest.out"
fi

# Section 3 - for real prediction: This is mainly imitated from braker pipeline. 
### ln -s ../../03.rnaseq_aln/02.map2unmsk/AllToCtg_hs2/hits.bam .
### ln -s ../../03.rnaseq_aln/02.map2unmsk/AllToCtg_hs2/hits.bam.bai .
##### make hints from BAM file
### ${dir_bin}/bam2hints  --intronsonly  --in=hits.bam  --out=bam2hints.temp.gff  2>bam2hints.0.err
##### add hints from BAM file to hints file
### cat bam2hints.temp.gff >> hintsfile.temp.gff
##### sort hints of type rnaseq
### cat hintsfile.temp.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 > hints.rnaseq.temp.sort.gff
##### join multiple hints
### perl $dir_aug/scripts/join_mult_hints.pl  < hints.rnaseq.temp.sort.gff > hintsfile.temp.gff 2> join_mult_hints.rnaseq.err
##### filter introns, find strand and change score to 'mult' entry
### perl $dir_braker/filterIntronsFindStrand.pl  input.genome.fa  hintsfile.temp.gff  --score  1>hintsfile.gff  2>filterIntronsFindStrand.err
##### rm hintsfile.temp.gff
##### splitting genome file in smaller parts for parallel execution of AUGUSTUS prediction
### perl $dir_aug/scripts/splitMfasta.pl  input.genome.fa  --outputpath=output_dir  --minsize=40000000  2> output_dir/splitMfasta.err
##### Run AUGUSTUS for each splitted genome fasta file like 'genome.split.*.fa'
### $dir_bin/augustus  --species=$orgName  --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH  --extrinsicCfgFile=$dir_braker/rnaseq.cfg  --alternatives-from-evidence=true  --hintsfile=hintsfile.gff  --UTR=off  --exonnames=on  --codingseq=on  output_dir/genome.split.1.fa  1>output_dir/augustus.1.gff  2>output_dir/augustus.1.stderr
##### .... Finish all predictions. 
##### Concatenating AUGUSTUS output files
### cat output_dir/augustus.1.gff >> output_dir/augustus.tmp.gff
### cat output_dir/augustus.2.gff >> output_dir/augustus.tmp.gff
##### .... Finish all files. 
##### Joining AUGUSTUS output
### perl $dir_aug/scripts/join_aug_pred.pl < output_dir/augustus.tmp.gff > output_dir/augustus.gff
##### Delete output_dir/augustus.*.gff files; 
##### Making a fasta file with protein sequences of output_dir/augustus.gff
### perl $dir_aug/scripts/getAnnoFasta.pl  output_dir/augustus.gff  --seqfile=input.genome.fa  2>output_dir/getAnnoFasta.augustus.err
##### Making a gtf file from output_dir/augustus.gff
### cat output_dir/augustus.gff | perl -ne 'if(m/\tAUGUSTUS\t/){print $_;}' | perl $dir_aug/scripts/gtf2gff.pl  --printExon --out=output_dir/augustus.gtf  2>output_dir/gtf2gff.augustus.gtf.err
##### Delete empty files by CMD: find  output_dir/  -empty
##### Done. 

cp -pr $AUGUSTUS_CONFIG_PATH/species/$orgName/ ./

