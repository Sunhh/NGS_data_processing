# ortho_tools — orthologs → alignments → gene trees → positive selection

From an all-vs-all protein comparison to 1-to-1 orthogroups, codon alignments, gene trees,
and a PAML branch-site test for positive selection. External tools: BLAST, an ortholog
clusterer (OrthoMCL / OrthoFinder), MUSCLE, ete3 (`ete3 build`), and PAML (`codeml`).

| script | role | step |
|--------|------|------|
| `mk_sep_blastp_shell.sh`, `replace_all_blast_file.sh` | split / assemble the all-vs-all blastp jobs | 0 |
| `filter_bp6_byTopScore.pl` | keep the top-scoring BLAST hits per query (`pep.bp6 > .slct`) | 0 |
| `01.ortho_list_from_orthoOut.pl` | OrthoMCL output → **1-to-1 orthogroups** (`.1to1_OGs`) | 1 |
| `ortho_cnt_c1_cmn.pl` | count OrthoFinder `Orthogroups.csv` groups by required taxa (`-required_taxaLis`) | 1 |
| `02.list_run_muscle.pl` | run **MUSCLE** per orthogroup (`-grp_list -pep_fas`) | 2 |
| `fmt_data_for_orthoAln.pl`, `sep_alnFas.pl`, `join_cdsFmt_faAln.pl` | prepare / split / join per-OG (CDS) alignments | 2 |
| `trim_faAln.pl` | trim a fasta alignment (`-in_fa -out`) | 2 |
| `cdsAln2bppAln.pl`, `cnvt_fa2nex.pl` | codon alignment → BPP / NEXUS format | 2 |
| `run_ete3_genetree.pl` | build gene trees with `ete3 build` (`-in_cog -in_prot -in_cds`, multi-threaded) | 3 |
| `combine_tree_for_treeAnnot.pl` | merge gene trees for annotation | 3 |
| `run_positive_selection.pl` | PAML **branch-site** test: null (`-in_tree0 -in_ctl0`) vs alternative (`-in_tree1 -in_ctl1`), LRT on a `#1`-labelled foreground branch | 4 |
| `cnvt_gffJnLoc_to_mcscanx.pl` | gff (joined-loc) → MCScanX gff (feeds the top-level synteny tools) | (synteny) |

---

## Workflow
```sh
# 0. all-vs-all blastp, then keep top hits, then cluster (OrthoMCL / OrthoFinder)
bash mk_sep_blastp_shell.sh ...           # split jobs
perl filter_bp6_byTopScore.pl pep.bp6 > pep.bp6.slct

# 1. 1-to-1 orthogroups
perl 01.ortho_list_from_orthoOut.pl all_orthomcl.out > all_orthomcl.out.1to1_OGs
#   (or count OrthoFinder groups by taxa: perl ortho_cnt_c1_cmn.pl -in_file Orthogroups.csv -required_taxaLis taxa)

# 2. align each orthogroup (MUSCLE) and trim
perl 02.list_run_muscle.pl -grp_list gene_grp_list -pep_fas all_pep.fas
perl trim_faAln.pl -in_fa OG_aln.fasta -out OG_aln.trim.fasta
#   codon-align helpers: join_cdsFmt_faAln.pl, cdsAln2bppAln.pl, cnvt_fa2nex.pl

# 3. gene trees
perl run_ete3_genetree.pl -in_cog all.ete3.cog -in_prot all.ete3.prot.fa -in_cds all.ete3.cds.fa -out_dir ete3_out/

# 4. positive selection (PAML branch-site LRT; mark the foreground branch with #1 in the tree)
perl run_positive_selection.pl -in_tree0 1by1.0.tree -in_tree1 1by1.1.tree -in_ctl0 1by1.0.ctl -in_ctl1 1by1.1.ctl ...
```
The 1-to-1 orthogroups also feed a species/SNP phylogeny (see `../../reseq_tools/phylo_tools/`
or `../ortho_tools` trees), and `cnvt_gffJnLoc_to_mcscanx.pl` feeds the synteny scripts in the
parent [`../`](../).

_Hand-written — `gen_script_index.sh` leaves it alone._
