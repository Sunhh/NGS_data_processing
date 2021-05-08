### Basic functions.
function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

soap_bin=/usr/local/bin/SOAPdenovo-127mer
graph_prefix='PG1_R02'
ctgK=91
mapK=91
cpuN=50
# multiK=127

# genomPara="-N 376191846"
genomPara=""

paraCtg="-M 2"

ctg_conf="PG1_ctg.cfg"
scf_conf="PG1_scf.cfg"

gapcloser_bin="/usr/local/bin/GapCloser"
gapcloser_lib="PG1_GC.cfg"
fa2close="${graph_prefix}.scafSeq"
fa2out="${graph_prefix}.GC"

# echo "[$(date)][Rec] Start $graph_prefix"
# $soap_bin all -s $scf_conf -o $graph_prefix -R -d 1 -K $ctgK -p $cpuN -a 200 -k $mapK -m $multiK -F -C 4 1>all.log.$graph_prefix 2>all.err.$graph_prefix

exe_cmd "$soap_bin pregraph -p $cpuN -s $ctg_conf -K $ctgK -R -d 1 -a 300 -o $graph_prefix 1>pregraph.log.$graph_prefix 2>pregraph.err.$graph_prefix"
exe_cmd "$soap_bin contig   -p $cpuN -s $ctg_conf -g $graph_prefix -R $paraCtg 1>contig.log.$graph_prefix 2>contig.err.$graph_prefix"
# echo "[$(date)][Msg] Please Edit file $graph_prefix.preGraphBasic file to fit mapK=$mapK"
exe_cmd "$soap_bin map      -p $cpuN -s $scf_conf -g $graph_prefix -k $mapK 1>map.log.$graph_prefix 2>map.err.$graph_prefix"
exe_cmd "$soap_bin scaff    -p $cpuN -g $graph_prefix -F $genomPara 1>scaff.log.$graph_prefix 2>scaff.err.$graph_prefix"

exe_cmd "$gapcloser_bin -a $fa2close -b $gapcloser_lib -o $fa2out -l 155 -p 31 -t $cpuN 1>gc.log.$graph_prefix 2>gc.err.$graph_prefix"

tsmsg "[Rec] All done."
