<!-- AUTOGEN gen_script_index.sh; do not edit — re-run the generator instead -->
# `project/watermelon_pan_phaseI` — script index

_24 Perl scripts. Auto-generated 2026-07-11; synopsis is from each script's own `perl $0` / `Usage: $0` line — run a script with `-h` or no args for full options._

| script | synopsis |
|---|---|
| `cnt_CDS_dup.pl` | cnt_CDS_dup.pl in.bn6 > in.bn6.ident_filteredNum |
| `cnt_gene_PAV_byDepCov_fromISbed.pl` | cnt_gene_PAV_byDepCov_fromISbed.pl 2_rdDep 90_coverage WM9_CalhounGray.CLpan.bam.gcovH.intersect > PAV.list |
| `cnt_gene_PAV_byDepCov.pl` | cnt_gene_PAV_byDepCov.pl 2_rdDep 90_coverage HS_CDS_loc.bed rdCov.bed > PAV.list |
| `cnt_TEIPRacc.pl` | cnt_TEIPRacc.pl potential_TE_IPRacc ipr_all6_tsv.TEprot.IPRacc.line > ipr_all6_tsv.TEprot.IPRacc.line.cnt |
| `cnvt_ext2nov_to_agp.pl` | cnvt_ext2nov_to_agp.pl CM.ext2nov.tbl > CM.ext2nov.agp |
| `comb_gene_PAV_byDepCov.pl` | comb_gene_PAV_byDepCov.pl CDScoveragePercent(0-100] file_list > summary.tbl |
| `convert_PAVtab_to_rmOlapTab.pl` | convert_PAVtab_to_rmOlapTab.pl final.wiRepre.fmt.PAV.1.rmOlap final.wiRepre.fmt.PAV > final.wiRepre.fmt.PAV.rmOlap |
| `dvd_grp_bySyn.pl` | dvd_grp_bySyn.pl -in_grp_mat ov2/comb.grp2.novl_loc -in_geneJnLoc input/all.gff3JnLoc -in_orthLis input/all.syn.tbl.orth -opref out |
| `extract_map_ratio_from_PEfstat_comb.pl` | extract_map_ratio_from_PEfstat_comb.pl comb.fstat_CLpan > comb.fstat_CLpan.tbl |
| `extract_N50.pl` | extract_N50.pl in.fa.N50 in_2.fa.N50 > out.table |
| `extract_novelOriFlankSeq.pl` | extract_novelOriFlankSeq.pl 100000(flank_length) map.dedup_rmcont_asm.tbl map.ID1_spec_procID_asmFn_asmPath out_prefix |
| `list_IPRacc.pl` | list_IPRacc.pl ipr_all6_tsv.TEprot.IPRacc > ipr_all6_tsv.TEprot.IPRacc.line |
| `map_accDedup_to_asm.pl` | map_accDedup_to_asm.pl Clean_UA.CD90.NR.fa.map2Src.tbl > map.dedup_rmcont_asm.tbl |
| `map_dedup_to_novel.pl` | map_dedup_to_novel.pl source.fa dedup.fa > dedup2source.tab |
| `pipe_for_functional_annotation.pl` | pipe_for_functional_annotation.pl annot.cfg |
| `ret_maker_abinit_gff3.pl` | ret_maker_abinit_gff3.pl ab_init.ID_list maker_all.gff3 > ab_init.match.gff3 |
| `ret_maker_abinit_gff3_simpleCDS.pl` | ret_maker_abinit_gff3_simpleCDS.pl ab_init.ID_list maker_all.gff3 > ab_init.match.gff3 |
| `rm_overlap_OGs.pl` | rm_overlap_OGs.pl final.wiDropped.cds.blk.gff3 final_rmQloc.wiRepre.fmt.PAV.5cols > final_rmQloc.wiRepre.fmt.PAV.rmOlap |
| `rm_Qloc_only_groups.pl` | rm_Qloc_only_groups.pl ov1/comb.grp2.novl_loc.syn.grp_tbl > ov1/comb.grp2.novl_loc.syn.grp_tbl.filtered |
| `rm_Qloc_ovl2pred.pl` | rm_Qloc_ovl2pred.pl -in_grp ov1/comb.grp2.novl_loc.syn.grp_tbl.filtered -in_blk ov1/all_type.cds.blk > ov1/comb.grp3 |
| `slct_OG_gene_pairs.pl` | slct_OG_gene_pairs.pl mcscan_run_p/ASM1002_97103_vs_ASM1003_PI537277.tbl.orth.gene_pairs > mcscan_run_p/ASM1002_97103_vs_ASM1003_PI537277.tbl.orth.gene_pairs.OG |
| `summary_covBp_byRdDep.pl` | summary_covBp_byRdDep.pl file_list > summary.tbl |
| `summary_gene_PAV_byDepCov.pl` | summary_gene_PAV_byDepCov.pl file_list > summary.tbl |
| `trim_gff_to_novel_ctg_region.pl` | trim_gff_to_novel_ctg_region.pl CL.ext2nov.tbl in_ext.gff3 > trimmed_ext.gff3 |
