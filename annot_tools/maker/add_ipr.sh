### Basic functions.
function exe_cmd {
  echo "[$(date)][CMD] $1"
  eval "$1"
  if [[ $? -eq "0" ]]
  then
    echo "[$(date)][CMD_done] $1"
  else
    echo "[$(date)][CMD_err] $1"
  fi
}

function tsmsg {
  echo "[$(date)][Msg] $1"
}

# mkPref="21TJS6_FalT1"
# gnPref="21TJS6"
mkPref=$1
gnFa=$2

IPRSCAN_HOME=/data/Sunhh/src/annotation/iprscan/interproscan-5.54-87.0/
PATH=/data/Sunhh/src/annotation/maker/maker/bin:$PATH

PL_dealTbl="deal_table.pl"
PL_retAbGff="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/ret_maker_abinit_gff3.pl"
PL_cvtAb2Mk="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/ret_makerGff_fromAbGff.pl"
PL_rmGffFa="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/rm_maker_fasta.pl"

exe_cmd "$IPRSCAN_HOME/interproscan.sh -iprlookup -goterms -pa -etra -cpu 24 -i ../${mkPref}.all.maker.non_overlapping_ab_initio.proteins.fasta -b m2 -f tsv,gff3,xml"
exe_cmd "$PL_dealTbl m2.tsv -kSrch_idx $IPRSCAN_HOME/../ipr.db_slct -kSrch_srcCol 3 | cut -f 1 | $PL_dealTbl -UniqColLine 0 > m2.tsv.wiIPR.mID"

# Get the gff3 file.
exe_cmd "$PL_retAbGff m2.tsv.wiIPR.mID ../3fal2_all.gff3 > m2.tsv.wiIPR.0.gff3"
### Convert ab to maker gff.
exe_cmd "$PL_cvtAb2Mk  ab  ../../../src_data/in_genomes/${gnFa}  m2.tsv.wiIPR.0.gff3"
# The final files are:
# gff3 : ab_maker_noFa.gff3
# protein sequence    : ab.all.maker.proteins.fasta
# transcript sequence : ab.all.maker.transcripts.fasta

### Combine gff files.
exe_cmd "$PL_rmGffFa  ../${mkPref}_maker.gff3 | grep -v ^# > comb_maker_noFa.gff3"
exe_cmd "cat ab_maker_noFa.gff3 >> comb_maker_noFa.gff3"

