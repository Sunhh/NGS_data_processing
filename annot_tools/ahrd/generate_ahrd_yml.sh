#!/usr/bin/bash

FN_protFa=$1
FN_out=$2
if [[ $FN_out == "" ]];
then
  echo "##########################################"
  echo "# Generate YML file for AHRD program."
  echo "bash $0  in.p.fa  ./v1/ahrd_output_v1.csv > ahrd_in_v1.yml"
  echo "##########################################"
  exit 1
fi
# ./v1/ahrd_output_v1.csv

#
# General parameters
#
FN_alnSprot="../01.blast2DB/o1.p2sprot.bp6"
FN_alnTremb="../01.blast2DB/o1.p2trembl.bp6"
FN_alnTair="../01.blast2DB/o1.p2tair10.bp6"

DIR_sampeList="/data/Sunhh/wmhifi/analysis/gene_functional_annotation/sample_list/"
FN_FaSprot="/data/Sunhh/database/db_fasta/uniprot/uniprot_sprot.fasta"
FN_FaTremb="/data/Sunhh/database/db_fasta/uniprot/uniprot_trembl.fasta"
FN_FaTair="/data/Sunhh/database/db_fasta/arabidopsis/Athaliana/TAIR10/others/TAIR10_pep_20110103_representative_gene_model_updated"


printf "proteins_fasta: $FN_protFa
token_score_bit_score_weight: 0.468
token_score_database_score_weight: 0.2098
token_score_overlap_score_weight: 0.3221
output: $FN_out

blast_dbs:
  swissprot:
    weight: 653
    description_score_bit_score_weight: 2.717061
    file: $FN_alnSprot
    database: $FN_FaSprot
    blacklist: $DIR_sampeList/blacklist_descline_sprot.txt
    filter: $DIR_sampeList/from_git/filter_descline_sprot.txt
    token_blacklist: $DIR_sampeList/blacklist_token.txt

  trembl:
    weight: 904
    description_score_bit_score_weight: 2.590211
    file: $FN_alnTremb
    database: $FN_FaTremb
    blacklist: $DIR_sampeList/blacklist_descline.txt
    filter: $DIR_sampeList/filter_descline_trembl.txt
    token_blacklist: $DIR_sampeList/filter_descline_trembl.txt

  tair:
    weight: 854
    description_score_bit_score_weight: 2.917405
    file: $FN_alnTair
    database: $FN_FaTair
"

echo '    fasta_header_regex: "^>(?<accession>[aA][tT][0-9mMcC][gG]\\d+(\\.\\d+)?)\\s+\\|[^\\|]+\\|\\s+(?<description>[^\\|]+)(\\s*\\|.*)?$"'
echo '    short_accession_regex: "^(?<shortAccession>.+)$"'
printf "    blacklist: $DIR_sampeList/blacklist_descline_tair.txt
    filter: $DIR_sampeList/filter_descline_tair.txt
    token_blacklist: $DIR_sampeList/blacklist_token.txt
"

