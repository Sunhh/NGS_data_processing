ksize=85
cpuN=10

function exe_cmd {
	echo "[$(date)][CMD] $1"
	eval $1
	echo "[$(date)][Rec] Done."
}

function tsmsg {
	echo "[$(date)]$1"
}

exe_quorum_mkdb='/data/share/src/Assemble/MaSuRCA/MaSuRCA-2.2.1/bin/quorum_create_database'
exe_jf2='/data/share/src/Assemble/MaSuRCA/MaSuRCA-2.2.1/bin/jellyfish-2.0'


in_dir="/data/Sunhh/apple_mtAsm/01.Clean/02.Q20/01.PE_Clean/Final"
in_seq=""
tsmsg "[Rec] All start."

#for inPref in C113_L1_pTr C113_L2_pTr
#do
#	in_seq="$in_dir/${inPref}_R1.fq $in_dir/${inPref}_R2.fq"
#	o_jf_pref="${inPref}_highQ_m${ksize}"
#	exe_cmd "$exe_jf2 count -m $ksize -s 8G -t $cpuN -F 4 -o ${o_jf_pref}.jf -c 8 -C $in_seq"
#	exe_cmd "$exe_jf2 histo ${o_jf_pref}.jf -h 9000000000 -t $cpuN -o ${o_jf_pref}.jf.histo"
#	exe_cmd "rm ${o_jf_pref}.jf"
#done

in_seq="$in_dir/C113_L1_pTr_R1.fq $in_dir/C113_L1_pTr_R2.fq"
in_seq="$in_seq $in_dir/C113_L2_pTr_R1.fq $in_dir/C113_L2_pTr_R2.fq"
o_jf_pref="C113_L1L2_pTr_highQ_m${ksize}"
exe_cmd "$exe_jf2 count -m $ksize -s 8G -t $cpuN -F 4 -o ${o_jf_pref}.jf -c 8 -C $in_seq"
exe_cmd "$exe_jf2 histo ${o_jf_pref}.jf -h 9000000000 -t $cpuN -o ${o_jf_pref}.jf.histo"
# exe_cmd "rm ${o_jf_pref}.jf"


tsmsg "[Rec] All done."


