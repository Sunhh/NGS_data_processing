<!-- AUTOGEN gen_script_index.sh; do not edit — re-run the generator instead -->
# `annot_tools` — script index

_36 Perl scripts. Auto-generated 2026-07-11; synopsis is from each script's own `perl $0` / `Usage: $0` line — run a script with `-h` or no args for full options._

| script | synopsis |
|---|---|
| `add_tag_to_gffID.pl` | add_tag_to_gffID.pl tag in.fmt.gff3 |
| `augustus.accuracy_calculator.pl` | augustus.accuracy_calculator.pl augusts_test1.stdout augusts_test2.stdout > all_accuracy.tbl |
| `clean_pasa_med_files.pl` | clean_pasa_med_files.pl in_transcript.fa pasa_dbName compreh_prefix |
| `cnvt_infernal_tbl.pl` | cnvt_infernal_tbl.pl P1Genom_V1.scf.fa.rRNA.tbl > P1Genom_V1.scf.fa.rRNA.tbl.tab |
| `cnvt_maker2aug_gff3.pl` | cnvt_maker2aug_gff3.pl in.maker_good.gff3 |
| `cnvt_spaln2makerAln_prot_gff3.pl` | cnvt_spaln2makerAln_prot_gff3.pl in.spaln.gff3 |
| `cnvt_uniprot_dat2fa.pl` | cnvt_uniprot_dat2fa.pl uniprot_sprot_plants.dat > uniprot_sprot_plants.fasta |
| `copy_species.pl` | newspecies.pl |
| `createAugustusJoblist.pl` | This is copied from augustus_scripts. |
| `find_complete_prot_byBlastp.pl` | find_complete_prot_byBlastp.pl -prot_qry incomplete_proteins.fasta -opref out_prefix |
| `fitGff_aug2maker.pl` | fitGff_aug2maker.pl aug1.out > aug1_withID_for_gff3_merge.gff |
| `fix_1bpLoc_by_zff2Gb.pl` | _(no synopsis found)_ |
| `get_gff_byScfID.pl` | get_gff_byScfID.pl -scfID scfID -suff |
| `get_sameCDSGff.pl` | get_sameCDSGff.pl -srcGff maker_r1.gff3 -idxGff maker_r2.gff3 > maker_r1_same2r2.gff3 |
| `good_model_from_gff3.pl` | good_model_from_gff3.pl -seedGff prev_good.gff3 |
| `intron2exex.pl` | make a fasta file with exon-exon sequences |
| `join_b2g_annot.pl` | join_b2g_annot.pl -gene_list gene_id_V1p7 b2g.final.annot > b2g.final.annot.perGene |
| `keep_nonRedundant_list.pl` | _(no synopsis found)_ |
| `merge_blast_xml.pl` | merge_blast_xml.pl r9GoodNameSub_00001.fasta.toNr.xml r9GoodNameSub_00002.fasta.toNr.xml |
| `pasa_gff_to_alnGff.pl` | pasa_gff_to_alnGff.pl -addTag |
| `pipe_get_complete_models.pl` | pipe_get_complete_models.pl out_prefix input.gff3 input.p.fas db_diamond |
| `predictByAug_rna2genome.pl` | predictByAug_rna2genome.pl -stepLis 123 [-continue] |
| `rename_by_GffJnLoc.pl` | rename_by_GffJnLoc.pl prefix_of_gene [-start_list in_list.chrN_nextGeneN] input.gff.JnLoc > map.oriMID_newMID_oriGID_newGID |
| `replace_blast_asn_db.pl` | replace_blast_asn_db.pl new_db_path input_raw.asn > new.asn |
| `rm_ExactDup_gene_model.pl` | rm_ExactDup_gene_model.pl dup_aug1_aug2.gff3 > nonDup_aug1_aug2.gff3 |
| `rm_ovlap_gene_model.pl` | rm_ovlap_gene_model.pl -srcGff gene_models.gff -idxGff repeat_loc.gff |
| `rmRedunt_inputNucl.pl` | rmRedunt_inputNucl.pl -nucl_qry redundant_nucleotides.fasta -opref out_prefix |
| `rmRedunt_inputProt.pl` | rmRedunt_inputProt.pl -prot_qry incomplete_proteins.fasta -opref out_prefix |
| `rmRepGff_withBadTarget.pl` | rmRepGff_withBadTarget.pl ID_class.not4Annot wm97pb_v2ID.scf.fa.mskAll_4maker.gff > wm97pb_v2ID.scf.fa.mskAll_4maker_clean.gff |
| `rmShrtExon_spaln_prot2genom.pl` | rmShrtExon_spaln_prot2genom.pl CM_CS_WM97_Arab_SprotPln.spaln_M4.fmt.gff3 |
| `run_spaln_prot2genom.pl` | run_spaln_prot2genom.pl -inFa protein.fa [-inFa in_prot_2.fa] -db db_genome |
| `set_stopCodonFreq.pl` | set_stopCodonFreq.pl etrain_rawGenes_findErr.std augustus/config/species/orgName/orgName_parameters.cfg |
| `simplify_gff3.pl` | simplify_gff3.pl in.gff3 > simple_IDs.gff3 |
| `slct_maker_gff3.pl` | slct_maker_gff3.pl mRNA_name.list maker.gff3 |
| `snap_good_wrn_by_valid.pl` | Read in genome.fat_valid file, and find good models. |
| `zff2augustus_gbk.pl` | 20150203 There is a bug when it meets a single base exon. It will set the strand to '-' without judging. |
