<!-- AUTOGEN gen_script_index.sh; do not edit — re-run the generator instead -->
# `annot_tools/liftoff_tools` — script index

_28 Perl scripts. Auto-generated 2026-07-11; synopsis is from each script's own `perl $0` / `Usage: $0` line — run a script with `-h` or no args for full options._

| script | synopsis |
|---|---|
| `blk2bed.pl` | blk2bed.pl output/CX.cds.blk.CL > output/CX.cds.blk.CL.bed |
| `blk2gff.pl` | blk2gff.pl output/CX.cds.blk.CL > output/CX.cds.blk.CL.gff3 |
| `chk_only_pan.pl` | chk_only_pan.pl comb.grp2.novl_loc.wiRepre.fmt > comb.grp2.novl_loc.wiRepre.fmt.ifOnlyPan |
| `cnt_R2Q_liftoff_info.pl` | cnt_R2Q_liftoff_info.pl R_own.gff3 Q_own.gff3 R_to_Q.liftoff.gff3 > R_to_Q.info.tbl |
| `cnvt_genemap_to_QlocBlk.pl` | cnvt_genemap_to_QlocBlk.pl genomic.fa output/mapCDS.CLpan.to.CMpan.tbl.gene_map > output/mapCDS.CLpan.to.CMpan.tbl.gene_map.Qcds.blk |
| `cnvt_genePAV_to_grpPAV.pl` | cnvt_genePAV_to_grpPAV.pl CL comb.grp2.novl_loc.tbl output/CL.mat.pav_rd1_cov50.good4pav > output/comb.grp2.CL_PAV |
| `cnvt_gff_to_cdsBed.pl` | cnvt_gff_to_cdsBed.pl out_prefix in.gff3 |
| `filter_R2Q_liftoff_tbl.pl` | _(no synopsis found)_ |
| `fit_gff_4igv.pl` | fit_gff_4igv.pl in.gff3 > in.rmGene.gff3 |
| `fmt_gff_trim2CDS.pl` | fmt_gff_trim2CDS.pl WCGv2.chr.gff3 > WCGv2.chr.trim2CDS.gff3 |
| `fmt_grp_by_spec.pl` | fmt_grp_by_spec.pl input/grp_tag comb.grp1 > comb.grp1.fmt |
| `grp2single.pl` | _(no synopsis found)_ |
| `info_bedtools_intersect.pl` | info_bedtools_intersect.pl bedtools_intersect_wao.out > bedtools_intersect_wao.out.tbl |
| `info_liftoffGff.pl` | info_liftoffGff.pl output/mapCDS.ASM1003_PI537277.to.ASM1002_97103.gff3 > output/mapCDS.ASM1003_PI537277.to.ASM1002_97103.gff3.tbl |
| `prepare_gff3_to_blk.pl` | prepare_gff3_to_blk.pl CLpan CLpan.trim2CDS.gff3.JnLoc > CLpan.trim2CDS.blk |
| `prepare_input.pl` | prepare_input.pl pref_list |
| `psl_to_geneLoc.pl` | psl_to_geneLoc.pl sp_tag in.psl > in.psl.gene_loc |
| `remove_ovl_loc.pl` | remove_ovl_loc.pl comb.grp2 > comb.grp2.novl_loc |
| `ret_good_model_from_liftoff_gff3.pl` | ret_good_model_from_liftoff_gff3.pl liftoff_out.gff3 > filtered.gff3 |
| `retrieve_liftoff_pairs.pl` | retrieve_liftoff_pairs.pl from_tag to_tag output/mapCDS.CCpan.to.CLpan.tbl |
| `retrieve_QlocSeq_fromBlk.pl` | retrieve_QlocSeq_fromBlk.pl genomic.fa output/Qcds.CX.blk.CL > output/Qcds.CX.blk.CL.fa |
| `retrieve_QlocSeq.pl` | retrieve_QlocSeq.pl genomic.fa output/mapCDS.CLpan.to.CMpan.tbl.gene_map > output/mapCDS.CLpan.to.CMpan.tbl.gene_map.Qcds.fa |
| `rm_head_partial_frame_inGff.pl` | rm_head_partial_frame_inGff.pl in.gff3 > in_fix.gff3 |
| `run_liftoff_FromTo_pan_para.pl` | run_liftoff_FromTo_pan_para.pl \\ |
| `run_liftoff_FromTo_pan.pl` | run_liftoff_FromTo_pan.pl from_XX to_XX |
| `run_liftoff_FromTo.pl` | run_liftoff_FromTo.pl from_XX to_XX |
| `simple_group_pairs.pl` | simple_group_pairs.pl output/mapCDS.C?pan.to.C?pan.tbl.gene_map > comb.grp1 |
| `slct_best_psl.pl` | slct_best_psl.pl in.psl > in.psl.gene_pair |
