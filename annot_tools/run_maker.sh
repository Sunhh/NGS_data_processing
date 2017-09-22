### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

cpuN=50
ctl_bopts="maker_bopts.ctl"
ctl_exe="maker_exe.ctl"

roundN=2
ctl_opts="maker_opts.r${roundN}.ctl"
exe_cmd "mpiexec -n $cpuN maker $ctl_opts $ctl_bopts $ctl_exe 1>maker.${roundN}-log 2>maker.${roundN}-err"

#roundN=3
#ctl_opts="maker_opts.r${roundN}.ctl"
#exe_cmd "mpiexec -n $cpuN maker $ctl_opts $ctl_bopts $ctl_exe 1>maker.${roundN}-log 2>maker.${roundN}-err"

# Get result. 
dir_maker='PG1All_v2_Scf.unmsk.maker.output'
db_idx='PG1All_v2_Scf.unmsk_master_datastore_index.log'
path_gffMerge=gff3_merge
path_fasMerge=fasta_merge

cd $dir_maker/
exe_cmd "mkdir ../result/r$roundN"
exe_cmd "$path_gffMerge -d $db_idx -g -o ../result/r${roundN}/r${roundN}_maker.gff3"
exe_cmd "$path_gffMerge -d $db_idx -n -o ../result/r${roundN}/r${roundN}_all.gff3"
exe_cmd "$path_fasMerge -d $db_idx -o ../result/r${roundN}/r${roundN}"
cd ../

# ipr_update_gff
# https://groups.google.com/forum/#!msg/maker-devel/VaoXWlGHOjs/kbh0YDl1b5gJ

# GlimmerHMM
# http://ccb.jhu.edu/software/glimmerhmm/man.shtml#spec_org

