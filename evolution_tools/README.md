# evolution_tools — comparative & evolutionary genomics

Synteny, orthology, structural variants, population structure, gene-family evolution and
selection between genomes/accessions. Each **analysis** lives in a subdirectory; the loose
top-level scripts cover **synteny plotting, RBH orthologs, and tree drawing**.

## Analysis pipelines (subdirectories)
| goal | where |
|------|-------|
| **Structural variants** from whole-genome alignment (AnchorWave/NucDiff → VCF) | [`SV_detection/`](SV_detection/) |
| **Copy-number variation** | [`copy_number_var/`](copy_number_var/) |
| **Orthology / gene trees / positive selection** (OrthoFinder → MUSCLE → ete3, PAML) | [`ortho_tools/`](ortho_tools/) |
| **Gene-family expansion / contraction** | [`expansion_tools/`](expansion_tools/) |
| **Population structure** (ADMIXTURE/STRUCTURE-style) | [`structure/`](structure/) |
| **Assembly-to-assembly comparison** | [`compare_assemblies/`](compare_assemblies/) |
| **MUMmer** helpers | [`mummer_tools/`](mummer_tools/) |
| **VCF-table** utilities | [`vcf_tab/`](vcf_tab/) |

---

## Top-level scripts

**Synteny**
| script | does |
|--------|------|
| `prepare_SynChro.pl` | prepare protein + gff input for SynChro (`prot.fa prot_chr.gff3 outPref badIDs`) |
| `follow_sibelia.pl` | post-process Sibelia synteny-block output (`-task ...`) |
| `cvt_mscGff_scf2chr.pl` | lift an MCScanX gff from scaffold to chromosome coordinates via an AGP (`-scf_gff -agp`) |
| `draw_syn_dotplot.pl` | synteny dot-plot between two scaffold sets (`-scfLis_x -scfLis_y -in_gff mcs.gff -in_aln mcs.collinearity`) |
| `plot_syn.pl` | ribbon/collinearity synteny plot between genomes |
| `plot_syn_bk.pl` | block-oriented synteny plot (a separate variant of `plot_syn.pl`) |

**Orthology (reciprocal best hits)**
| script | does |
|--------|------|
| `rbh_byBp6.pl` | reciprocal best hits from BLAST tabular (`-in_bp6`; filter `-min_similarity -max_lenDiffR`) |
| `rbh_inBlock.pl` | RBH restricted to MCScanX collinearity blocks (wraps `rbh_byBp6.pl`) |

**Tree plotting** (`Rscript`, needs **ggtree**)
| script | does |
|--------|------|
| `plot_color_tree.r` | draw a tree with branch/tip colouring |
| `plot_label_color_tree.r` | as above, with tip labels |

_Hand-written — `gen_script_index.sh` leaves it alone. Subdir `README.md` files hold the
per-analysis details; `../SCRIPTS_INDEX.md` is the flat repo-wide index._
