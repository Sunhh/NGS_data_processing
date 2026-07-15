# expansion_tools — gene-family expansion / contraction (CAFE)

Test which gene families have significantly expanded or contracted along a species tree with
**CAFE**, starting from orthogroups (see [`../ortho_tools/`](../ortho_tools/)) and an
ultrametric species tree, then extract and annotate the significant families.

External tool: **CAFE** (Computational Analysis of gene Family Evolution).

| script | role | step |
|--------|------|------|
| `01.prepare_ortho_to_tbl.pl` | OrthoMCL output → a plain family × gene table | 1 |
| `01.prepare_cafe_tab.pl` | OrthoMCL output → CAFE count table (family × taxon gene counts; `-taxa_list`) | 1 |
| `01.clean_nwk.pl` | tidy a Newick species tree for CAFE (ultrametric input) | 1 |
| `03.replace_geneID_in_orthomcl.pl` | rename gene IDs in the orthomcl output (`-in_name_list -in_orthomcl`) | 1 (optional) |
| `02.cafe_to_grp.pl` | CAFE report → per-node expanded/contracted family groups (`-cafe_report`) | 3 |
| `04.get_expansion_tab.pl` | CAFE solution table → expansion table (`-ratio -bigCol -smallCol`) | 3 |
| `05.add_desciption_to_OGcsv.pl` | add functional descriptions (AHRD) to `Orthogroups.csv` (`-tab_desc`) | 4 |
| `jn_gene_byIPR.pl` | group genes by InterPro (IPR) annotation from an iprscan tsv | 4 |

---

## Workflow
```sh
# 1. CAFE input: gene-count table + ultrametric species tree
perl 01.prepare_cafe_tab.pl all_orthomcl.out -taxa_list taxa_list > cafe_input.tab
perl 01.clean_nwk.pl species.nwk > species.ultrametric.nwk

# 2. run CAFE (external) on cafe_input.tab + the tree  ->  input.cafe / *.report

# 3. extract significant expansions / contractions
perl 02.cafe_to_grp.pl -cafe_report input.cafe
perl 04.get_expansion_tab.pl cafe_input_grp01.report.cafe.sol.tab -ratio 1 -bigCol 3 -smallCol 4 \
  > cafe_input_grp01.report.cafe.sol.tab.sol_big

# 4. annotate the families
perl 05.add_desciption_to_OGcsv.pl -tab_desc P1_ahrd.2c Orthogroups.csv > Orthogroups.csv.desc
perl jn_gene_byIPR.pl in.iprV5.tsv
```
Orthogroups come from `../ortho_tools/` (OrthoMCL / OrthoFinder).

_Hand-written — `gen_script_index.sh` leaves it alone. (`05.add_desciption_*` keeps its
original filename spelling to avoid breaking existing references.)_
