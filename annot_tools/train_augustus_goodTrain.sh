### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dir_aug='/Data/Sunhh/src/annotation/maker/maker/exe/augustus/'
dir_cfg="${dir_aug}/config/"
dir_bin="${dir_aug}/bin/"

dir_braker='/Data/Sunhh/src/annotation/braker/BRAKER_v2.0.4/'

export PATH="$dir_aug/bin/:$PATH"


pl_zff2augustus=$HOME/tools/github/NGS_data_processing/annot_tools/zff2augustus_gbk.pl
pl_fixZff2Aug=$HOME/tools/github/NGS_data_processing/annot_tools/fix_1bpLoc_by_zff2Gb.pl
pl_setStopCod=$HOME/tools/github/NGS_data_processing/annot_tools/set_stopCodonFreq.pl
pl_accuracy=$HOME/tools/github/NGS_data_processing/annot_tools/augustus.accuracy_calculator.pl

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
# ln -s ../../../07.predBySnap/01.train_snap/03.train_from_r4finalUse/export.dna .
# ln -s ../../../07.predBySnap/01.train_snap/03.train_from_r4finalUse/export.ann .
orgRef='RoseR1Pasa'
# orgRef='wm97pbv0Brk0'
orgName='RoseMakerR4finalUse'
testNum=1000
# in_genom_fa=export.dna
# in_genom_zff=export.ann
# in_genom_pep=export.aa
# in_raw_gff=r5_maker_good_noFa.gff3
# in_genom_fa=P1All.scf.fa
# cp -p /Data/Sunhh/seldom/forWushan/rose_genome/gene_prediction/03.rnaseq_aln/04.PASA_pipe/01.denovo_and_GG/filter_by_ProtDB/cand_merged_finalUse.genomeID.gb .
in_raw_gb="genes_for_${orgName}.gb"
in_gb=train.gb

fn_orgParaCfg=$dir_aug/config/species/$orgName/${orgName}_parameters.cfg

### Generate list to remove redundancy by proteins similarity
#exe_cmd "perl $pl_rmRedFromProt -minIdentity 70 -prot_qry $in_genom_pep -opref ${in_genom_pep}.rmR"
#para_rmLis="-rmLis=${in_genom_pep}.rmR.bp6.redund_list"
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
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${dir_cfg} $in_raw_gb.fix 1>etrain_findErr.std 2>etrain_findErr.err"
exe_cmd "cat etrain_findErr.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb.fix > $in_gb"
exe_cmd "perl $pl_randSplit $in_gb $testNum"
exe_cmd "perl $pl_setStopCod etrain_findErr.std $fn_orgParaCfg"
##echo "vim /workdir/laopopo/src/Annot/maker/maker/exe/augustus/config/species/$orgName/${orgName}_parameters.cfg"
##echo "Set stopCodonExcludedFromCDS to false , because our self genbank includes stop codon."
#exit; 

# Section 2
### Training without optimize. If the new species_model comes from a trained species, this step is not requried. 
###   And if the training dataset is bad, say redundant or with some other problem, the further training will make the model worse!!! 
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${dir_cfg} ${in_gb}.train 1>etrain_trainRaw.std 2>etrain_trainRaw.err"
# Test the first model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName  --AUGUSTUS_CONFIG_PATH=${dir_cfg}  ${in_gb}.test > goodTrain_dbGood_firsttest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttest.out"
exe_cmd "perl $pl_accuracy goodTrain_dbGood_firsttest.out"
### Till now, I have never successfully trained UTR model, so I give up. 
### Training without optimize. UTR=on
#exe_cmd "$exe_etrain --species=$orgName --UTR=on --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainRawUtr.std 2>etrain_trainRawUtr.err"
# Test the first model by an ab initial prediction
#exe_cmd "$exe_augustus --species=$orgName --UTR=on ${in_gb}.test > goodTrain_dbGood_firsttestUtr.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttestUtr.out"

### Training with optimize. 
# exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
# exe_cmd "perl $pl_optPara --cpus=50 --kfold=25 --UTR=on --species=$orgName  ${in_gb}.train 1>opt_train_utr.std 2>opt_train_utr.err"
exe_cmd "perl $pl_optPara --cpus=100 --rounds=5 --AUGUSTUS_CONFIG_PATH=${dir_cfg} --species=$orgName  --onlytrain=${in_gb}.train  ${in_gb}.test  1>opt_train.std 2>opt_train.err"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false --AUGUSTUS_CONFIG_PATH=${dir_cfg} ${in_gb}.train 1>etrain_trainOpt.std 2>etrain_trainOpt.err"
## Test the second model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName  --AUGUSTUS_CONFIG_PATH=${dir_cfg}  ${in_gb}.test > goodTrain_dbGood_secondtest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_secondtest.out"
exe_cmd "perl $pl_accuracy goodTrain_dbGood_secondtest.out"


#exe_cmd "$exe_augustus --species=$orgName cegma_GoldDb.gff > goodTrain_dbCegma_firsttest.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbCegma_firsttest.out"

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
### $dir_bin/augustus  --species=$orgName  --AUGUSTUS_CONFIG_PATH=$dir_cfg  --extrinsicCfgFile=$dir_braker/rnaseq.cfg  --alternatives-from-evidence=true  --hintsfile=hintsfile.gff  --UTR=off  --exonnames=on  --codingseq=on  output_dir/genome.split.1.fa  1>output_dir/augustus.1.gff  2>output_dir/augustus.1.stderr
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


