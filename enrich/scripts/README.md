<!-- AUTOGEN gen_script_index.sh; do not edit — re-run the generator instead -->
# `enrich/scripts` — script index

_13 Perl scripts. Auto-generated 2026-07-11; synopsis is from each script's own `perl $0` / `Usage: $0` line — run a script with `-h` or no args for full options._

| script | synopsis |
|---|---|
| `cnvt_GOobo_to_tab.pl` | cnvt_GOobo_to_tab.pl > gene_ontology_edit.obo.2018-05-01.tab |
| `enrich_IPR.pl` | enrich_IPR.pl background_IPR subset_geneList [padj_cutoff_0.05] |
| `enrich_keggPWY.pl` | enrich_keggPWY.pl background_pwy subset_geneList [padj_cutoff_0.05] |
| `enrichment_mine_fit.pl` | enrichment_mine_fit.pl background_IPR subset_list |
| `extend_GOannot_for_GOenrich.pl` | extend_GOannot_for_GOenrich.pl in_go.obo.tab in_GO.annot out_file |
| `extend_IPRannot_for_IPRenrich.pl` | extend_IPRannot_for_IPRenrich.pl pub_IPR-entry.list in_gene2IPR.combined out_file |
| `hs_enrich_old.pl` | hs_enrich_old.pl -inSubList subset_geneList -bgGOtab wm97pbV2ID_feiID_clean.GO_bg.tab -oboTab gene_ontology_edit.obo.2018-05-01.tab |
| `hs_enrich.pl` | hs_enrich.pl -inSubList subset_geneList -bgGOtab wm97pbV2ID_feiID_clean.GO_bg.tab -oboTab gene_ontology_edit.obo.2018-05-01.tab_info |
| `mk_bg_fromKeggPWY.pl` | mk_bg_fromKeggPWY.pl KEGG_PWY.txt -colN 4 -bg_geneList Cma_geneID_list > Cma_geneID_list.kgPWY_bg |
| `prepare_enrich.pl` | prepare_enrich.pl -bg_geneList genome_geneID_list -in keggPWYByKO.tab |
| `split_jnGO.pl` | split_jnGO.pl -in wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col |
| `stat_goslim.pl` | stat_goslim.pl -in_idlist in_geneID.list -opref out_prefix -ref_goslim goslim_plant.obo.20181129 -ref_bgobo gene_ontology_edit.obo.2018-05-01 -ref_bggaf ngsV1.GAFngsV1.GAF |
| `test_enrichGO.pl` | test_enrichGO.pl background_pwy subset_geneList [padj_cutoff_0.05] |
