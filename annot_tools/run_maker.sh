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

roundN=1
ctl_opts="maker_opts.r${roundN}.ctl"
exe_cmd "mpiexec -n $cpuN maker $ctl_opts $ctl_bopts $ctl_exe 1>maker.${roundN}-log 2>maker.${roundN}-err"

#roundN=2
#ctl_opts="maker_opts.r${roundN}.ctl"
#exe_cmd "mpiexec -n $cpuN maker $ctl_opts $ctl_bopts $ctl_exe 1>maker.${roundN}-log 2>maker.${roundN}-err"

# ipr_update_gff
# https://groups.google.com/forum/#!msg/maker-devel/VaoXWlGHOjs/kbh0YDl1b5gJ

# GlimmerHMM
# http://ccb.jhu.edu/software/glimmerhmm/man.shtml#spec_org

