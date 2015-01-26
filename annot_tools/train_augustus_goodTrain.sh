### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

dir_aug='/workdir/laopopo/src/Annot/maker/maker/exe/augustus'

export PATH="$PATH:$dir_aug/bin/"


pl_zff2augustus=$HOME/tools/github/NGS_data_processing/annot_tools/zff2augustus_gbk.pl
pl_gff2gb=$dir_aug/scripts/gff2gbSmallDNA.pl
pl_randSplit=$dir_aug/scripts/randomSplit.pl
pl_newSpec=$dir_aug/scripts/new_species.pl
pl_filtGene=$dir_aug/scripts/filterGenes.pl
pl_optPara=$dir_aug/scripts/optimize_augustus.pl

exe_etrain=etraining
exe_augustus=augustus

# Need export.dna and export.ann files. 
orgName='SPGr1Good_AED0'
testNum=200
#in_genom_fa=export.dna
#in_raw_gff=export_trimStop.gff3
in_raw_gb=export_wiStop.gb
in_gb=train.gb


# perl $pl_gff2gb $in_raw_gff $in_genom_fa 2000 $in_raw_gb
# exe_cmd "$exe_etrain --species=generic --stopCodonExcludedFromCDS=true $in_raw_gb 1>etrain_findStop.std 2>etrain_findStop.err"
# exe_cmd "cat etrain_findStop.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
# exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb > $in_gb"

exe_cmd "perl $pl_zff2augustus > $in_raw_gb"
exe_cmd "$exe_etrain --species=generic --stopCodonExcludedFromCDS=false $in_raw_gb 1>etrain_findErr.std 2>etrain_findErr.err"
exe_cmd "cat etrain_findErr.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > dropgenes.lst"
exe_cmd "perl $pl_filtGene dropgenes.lst $in_raw_gb > $in_gb"
exe_cmd "perl $pl_randSplit $in_gb $testNum"
exe_cmd "perl $pl_newSpec --species=$orgName"
#echo "vim /workdir/laopopo/src/Annot/maker/maker/exe/augustus/config/species/$orgName/${orgName}_parameters.cfg"
#echo "Set stopCodonExcludedFromCDS to false , because our self genbank includes stop codon."
#exe_cmd "perl $pl_optPara --cpus=50 --species=$orgName ${in_gb}.train 1>opt_train.std 2>opt_train.err"
exe_cmd "$exe_etrain --species=$orgName --stopCodonExcludedFromCDS=false ${in_gb}.train 1>etrain_train.std 2>etrain_train.err"

## Test the first model by an ab initial prediction
exe_cmd "$exe_augustus --species=$orgName ${in_gb}.test > goodTrain_dbGood_firsttest.out"
exe_cmd "grep -A 22 Evaluation goodTrain_dbGood_firsttest.out"

#exe_cmd "$exe_augustus --species=$orgName cegma_GoldDb.gff > goodTrain_dbCegma_firsttest.out"
#exe_cmd "grep -A 22 Evaluation goodTrain_dbCegma_firsttest.out"

