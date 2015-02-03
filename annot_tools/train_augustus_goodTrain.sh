### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval "$1"
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dir_aug='/data/Sunhh/src/Annot/maker/maker/exe/augustus'

export PATH="$PATH:$dir_aug/bin/"


pl_zff2augustus=$HOME/tools/github/NGS_data_processing/annot_tools/zff2augustus_gbk.pl
para_rmLis=""
pl_gff2gb=$dir_aug/scripts/gff2gbSmallDNA.pl
pl_randSplit=$dir_aug/scripts/randomSplit.pl
pl_newSpec=$dir_aug/scripts/new_species.pl
pl_filtGene=$dir_aug/scripts/filterGenes.pl
pl_optPara=$dir_aug/scripts/optimize_augustus.pl

pl_rmRed=$HOME/tools/github/NGS_data_processing/annot_tools/keep_nonRedundant_list.pl
pl_dealTbl=$HOME/tools/github/NGS_data_processing/deal_table.pl


exe_etrain=etraining
exe_augustus=augustus

# Need export.dna and export.ann files. 
orgName='SPGr2Good_AED0'
testNum=200
in_genom_fa=export.dna
in_genom_zff=export.ann
in_genom_pep=export.aa
#in_raw_gff=export_trimStop.gff3
in_raw_gb=export_wiStop.gb
in_gb=train.gb


### Generate list to remove redundancy by proteins similarity
exe_cmd "makeblastdb -dbtype prot -in $in_genom_pep"
exe_cmd "blastp -task blastp -query $in_genom_pep -db $in_genom_pep -num_threads 10 -out ${in_genom_pep}.self.bp6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand' "
self_bp6="${in_genom_pep}.self.bp6"
self_bp6_slct="$self_bp6.slct"
exe_cmd "awk ' \$1 != \$2 && \$3 >= 70 && \$4 >= \$13 * 0.7 ' ${in_genom_pep}.self.bp6 > ${in_genom_pep}.self.bp6.slct "
exe_cmd "perl $pl_rmRed ${in_genom_pep}.self.bp6.slct > ${in_genom_pep}.self.bp6.redund_list"
para_rmLis="-rmLis=${in_genom_pep}.self.bp6.redund_list"
exe_cmd "rm ${in_genom_pep}.phr ${in_genom_pep}.pin ${in_genom_pep}.psq"


### Prepare for augustus and etrain. 
# perl $pl_gff2gb $in_raw_gff $in_genom_fa 2000 $in_raw_gb
# exe_cmd "$exe_etrain --species=generic --stopCodonExcludedFromCDS=true $in_raw_gb 1>etrain_findStop.std 2>etrain_findStop.err"
# exe_cmd "cat etrain_findStop.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
# exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb > $in_gb"

exe_cmd "perl $pl_zff2augustus $para_rmLis > $in_raw_gb"
exe_cmd "$exe_etrain --species=generic --stopCodonExcludedFromCDS=false $in_raw_gb 1>etrain_findErr.std 2>etrain_findErr.err"
exe_cmd "cat etrain_findErr.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb > $in_gb"
exe_cmd "perl $pl_randSplit $in_gb $testNum"
exe_cmd "perl $pl_newSpec --species=$orgName"
#echo "vim /workdir/laopopo/src/Annot/maker/maker/exe/augustus/config/species/$orgName/${orgName}_parameters.cfg"
#echo "Set stopCodonExcludedFromCDS to false , because our self genbank includes stop codon."

### Training without optimize. 
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainRaw.std 2>etrain_trainRaw.err"
## Test the first model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName ${in_gb}.test > goodTrain_dbGood_firsttest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttest.out"

### Training with optimize. 
exe_cmd "perl $pl_optPara --cpus=50 --species=$orgName ${in_gb}.train 1>opt_train.std 2>opt_train.err"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_trainOpt.std 2>etrain_trainOpt.err"
## Test the second model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName ${in_gb}.test > goodTrain_dbGood_secondtest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_secondtest.out"


#exe_cmd "$exe_augustus --species=$orgName cegma_GoldDb.gff > goodTrain_dbCegma_firsttest.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbCegma_firsttest.out"

