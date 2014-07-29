#!/bin/bash

## Basic information.
args=("${@}")
myBasename=${0##*/}
path2Me=${0}

# Configuration

## input : user defined.
refFa="${args[0]}"
qryFa="${args[1]}"
dbTag="${args[2]}"
oPref="${args[3]}"
scfPref="${args[4]}"
# mkdir db/
if [[ -n "$refFa" && -n "$qryFa" && -f "$refFa" && -f "$qryFa" ]]
then
	echo "[$(date)][Rec] refFa=$refFa qryFa=$qryFa"
else
	echo "[$(date)][Err] Please check the input for refFa and qryFa."
	echo "[$(date)][CMD] bash $path2Me/$myBasename refFa qryFa dbTag oPref scfPref"
	exit 1
fi

if [[ -n "$dbTag" ]]
then
	echo "[$(date)][Rec] dbTag=$dbTag"
else
	dbTag="dbTag"
	echo "[$(date)][Rec] dbTag set to $dbTag"
fi

if [[ -n "$oPref" ]]
then
	echo "[$(date)][Rec] oPref=$oPref"
else
	oPref="oPref"
	echo "[$(date)][Rec] oPref set to $oPref"
fi

if [[ -n $scfPref ]]
then
	echo "[$(date)][Rec] scfPref=$scfPref"
else
	scfPref=$oPref
	echo "[$(date)][Rec] scfPref set to $scfPref"
fi

### Medium files. 
file_refClip="${refFa}.clip"
file_qryClip="${qryFa}.clip"
file_refCtg="${oPref}_Ref.fitMugsy.fa"
file_qryCtg="${oPref}_Qry.fitMugsy.fa"
file_ctgFa="${oPref}_all.fitMugsy.fa"
file_rawMaf="${oPref}.maf"
file_fitMaf="${oPref}.fitMugsy.maf"
file_xmfa="${oPref}.fitMugsy.maf.xmfa"
pref_mugO="${oPref}_syn"
file_link="${pref_mugO}.paired.maf.link"
file_tbl="${oPref}.tbl"
file_scaf="${oPref}.scaf"
file_scafSeq="${oPref}.scafSeq"
file_drop="${file_scafSeq}.dropped"


### Scripts and softwares. 
pl_clipScaf="$HOME/tools/github/NGS_data_processing/assemble_tools/clip_scaf_end.pl"
pl_addTag2fa="$HOME/tools/github/NGS_data_processing/assemble_tools/add_tag_to_fsa.pl"
pl_formatMAF="$HOME/tools/github/NGS_data_processing/assemble_tools/format_maf_forMugsy.pl"
pl_maf2fasta="$HOME/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl"

sh_mugsyenv="$HOME/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh"
exe_mugsyWGA="$HOME/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA"
pl_pairedMAF="$HOME/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl"

pl_goodLink="$HOME/tools/github/NGS_data_processing/assemble_tools/good_link_fromMaf.pl"
pl_joinScaf="$HOME/tools/github/NGS_data_processing/assemble_tools/join_link.pl"

exe_lastdb="$HOME/bin/lastdb"
exe_lastal="$HOME/bin/lastal"
exe_paraFa="$HOME/bin/parallel-fasta"
exe_lastSp="$HOME/bin/last-split"

echo "[$(date)] Clip input scaffolds."
perl $pl_clipScaf $refFa > $file_refClip
perl $pl_clipScaf $qryFa > $file_qryClip

echo "[$(date)] Format sequences."
perl $pl_addTag2fa Ref $file_refClip > $file_refCtg
perl $pl_addTag2fa Qry $file_qryClip > $file_qryCtg
cat $file_refCtg $file_qryCtg > $file_ctgFa

echo "[$(date)] lastdb"
$exe_lastdb -w 5 -m111111111111111111111111111110 $dbTag $file_refCtg 

echo "[$(date)] lastal"
$exe_paraFa "$exe_lastal -q3 -e30 $dbTag | $exe_lastSp " < $file_qryCtg > $file_rawMaf

echo "[$(date)] Format alignments."
perl $pl_formatMAF -specs "," $file_rawMaf -out $file_fitMaf
perl $pl_maf2fasta $file_fitMaf > $file_xmfa

echo "[$(date)] Run Mugsy"
source $sh_mugsyenv                                                                                                                                                     $exe_mugsyWGA --outfile $pref_mugO --seq $file_ctgFa --aln $file_xmfa --distance 2000 --minlength 50 1>stdout.${pref_mugO} 2>stderr.${pref_mugO}
perl $pl_pairedMAF ${pref_mugO}.maf > ${pref_mugO}.paired.maf

echo "[$(date)] Looking for links"
perl $pl_goodLink ${pref_mugO}.paired.maf -minClustLen 2000 1>${file_link} 2>${file_link}.err

echo "[$(date)] Build scaffolds"
perl $pl_joinScaf ${file_link} -inCtgFa $file_ctgFa -outLnkTbl ${file_tbl} -outScfLnk ${file_scaf} -outScfFas ${file_scafSeq} -scfPref ${scfPref} -dropScfIn ${file_drop}

# echo "[$(date)] Get"


echo "[$(date)] All done."


