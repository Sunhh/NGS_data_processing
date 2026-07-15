# coexp — WGCNA co-expression modules & trait association

Build a WGCNA co-expression network from an expression matrix, then associate the modules
with sample traits. The network is built **once** (`run_wgcna.r`) and can be associated
with several trait sets (`down_phenoAssoc.r`) without rebuilding.

Needs `R` with **WGCNA** (and `pheatmap` + `RColorBrewer` for the module-comparison heatmap).

## Method notes (WGCNA is still standard — a few current-practice defaults)
WGCNA remains a standard method for bulk RNA-seq co-expression; its core algorithm
(soft-thresholding → TOM → dynamic tree cut → module eigengenes) is unchanged. These
scripts default to current best practice; override if you need the older behaviour:
- **`--transform log2`** (default): build the network on `log2(expr+0.01)`, **not** raw
  FPKM/RPKM — WGCNA authors advise log/variance-stabilised input. `--transform none`
  reproduces the old raw-value behaviour.
- **`--corType bicor`** (default): biweight midcorrelation is more robust to outliers than
  Pearson. `--corType pearson` for the classic correlation.
- **`--networkType signed`** (default): signed networks separate positively/negatively
  correlated genes; generally preferred over `unsigned`.
- **Sample size**: correlation networks want enough observations — aim for ≳15–20
  samples/conditions. (Group-mean expression reduces the effective N.)

For heavily automated module analysis + enrichment, `CEMiTool` is a modern WGCNA-based
alternative; for single-cell/spatial data see `hdWGCNA`. For bulk data with full control,
classic WGCNA (here) is appropriate.

## Scripts
| script | role |
|--------|------|
| `run_wgcna.r` | build the network → module assignment, eigengenes (`*_network_ForPhenoAssoc.RData`), per-gene KME table, diagnostic pdfs, TOM + Cytoscape export |
| `down_phenoAssoc.r` | module ↔ trait association from the network RData → labelled heatmap (pdf) + `cor`/`pval` table (tsv) |
| `dist_of_twoKME.pl` + `heatmap_by_mod_dist.r` | compare two KME tables (e.g. signed vs unsigned) → similarity matrix → distance heatmap |
| `add_abs_toKME.pl` | append an `absKME` column to a bare KME table (the `run_wgcna.r` KME already has one) |

---

## (1) Inputs
- **Expression** (`--expr`): row1 header, col1 = gene ID, other columns = per-sample
  expression (a mean-per-group matrix is fine); missing values as `-`.
- **Traits** (`--pheno`): row1 header, col1 = sample ID (matching the expression columns),
  other columns = traits (continuous or binary 0/1; `NA` allowed).

## (2) Build the network
```sh
Rscript run_wgcna.r --expr input/dat1_rpkmMean --outdir wgcna_dat1 --prefix dat1
# key options: --networkType signed|unsigned  --corType bicor|pearson
#              --transform log2|none  --softPower <int>  --minModuleSize 30  --mergeCutHeight 0.25
```
Outputs under `wgcna_dat1/`: `dat1_network_ForPhenoAssoc.RData` (modules + eigengenes for
step 3), `dat1_KME.txt` (per-gene module membership; sorted by module then |KME|, with an
`absKME` column), `dat1_sftThreshold.pdf`, `dat1_dendroColors.pdf`, `dat1_adjHeatmap.pdf`,
`dat1_ME_barplot/`, `dat1_cytoscape/` (edge/node files per module), and the TOM RData.

## (3) Associate modules with traits
```sh
Rscript down_phenoAssoc.r --network wgcna_dat1/dat1_network_ForPhenoAssoc.RData \
  --pheno input/dat1_pheno --out wgcna_dat1/dat1 --corType bicor
# add --binaryTrait if the traits are 0/1 (sets bicor robustY=FALSE)
```
Writes `<out>.module_trait.heatmap.pdf` (module-trait + module-module heatmaps, `Sig`
flags at `--traitPval`, default 0.01) and `<out>.module_trait.tsv` (module × trait
correlation and p-value). Rerun with different `--pheno` / `--corType` on the same
network — no rebuild.

## (4) Optional — compare two networks / export to Cytoscape
Compare the module membership of, say, a signed vs an unsigned build:
```sh
perl dist_of_twoKME.pl signed wgcna_signed/dat1_KME.txt unsigned wgcna_unsigned/dat1_KME.txt > KME_sim.tab
Rscript heatmap_by_mod_dist.r KME_sim.tab KME_sim.heatmap.pdf
```
`wgcna_dat1/dat1_cytoscape/input-{edges,nodes}-<module>.txt` load directly into Cytoscape.

---

_See `run_example.sh` for a runnable demo on `input/dat1_*` (subsets to the top variable
genes so it finishes quickly). Hand-written — `gen_script_index.sh` leaves it alone._
