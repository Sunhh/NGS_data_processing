# phylo_tools — SNP phylogeny (random-sampled datasets → NJ tree)

Build a whole-genome-SNP phylogeny: draw a random SNP dataset, run a Neighbor-Joining tree
in **MEGA**, and relabel the tips. `generate_dataset.pl` automates the SNP-table → alignment
steps; the `.mao` is the MEGA analysis template; `replace_ID_in_nwk.pl` fixes tip names.

Needs `MEGA` (megacc) for the tree; `generate_dataset.pl` calls repo helpers
(`../rand_site_wiWind.pl`, `../cnvt_tools/cols2fas.pl`, `deal_table.pl`, `deal_fasta.pl`).

| file | role |
|------|------|
| `generate_dataset.pl` | SNP table → a random-site, per-taxa FASTA alignment for one replicate (`-dirID -setID` name the output dir) |
| `infer_NJ_nucleotide_bt500.mao` | MEGA analysis options: Neighbor-Joining, nucleotide, 500 bootstraps |
| `replace_ID_in_nwk.pl` | rewrite tip IDs in a Newick tree (`old2new_ID in.nwk > new.nwk`) |

---

## Workflow
1. Build a replicate alignment (random SNP sites within windows, subset to the taxa list,
   converted to FASTA and relabelled):
```sh
perl generate_dataset.pl -snp_tbl Acc131_mask.snp -opref Acc131_mask \
  -tax_list GrpList/grp52_list -dirID 01 -setID 01 -wLen 10000 -wNum 30
# -> 01_rand_01/Acc131_mask_set01.*.fa   (repeat with new -setID for more replicates)
```
2. Infer the NJ tree with MEGA-CC using the bundled options:
```sh
megacc -a infer_NJ_nucleotide_bt500.mao -d 01_rand_01/Acc131_mask_set01.use.snp.fa -o tree01
```
3. Relabel tips (e.g. accession IDs → display names):
```sh
perl replace_ID_in_nwk.pl old2new_ID tree01.nwk > tree01.named.nwk
```

## Method note
Neighbor-Joining (MEGA) on random SNP subsets is fast and fine for a quick overview, but for
a publication-grade SNP phylogeny prefer a **maximum-likelihood** tree with ascertainment-bias
correction — e.g. **IQ-TREE** (`-m GTR+ASC`) or **RAxML** — on the concatenated SNP alignment,
with ultrafast bootstrap. (An ML/IQ-TREE workflow is used elsewhere in this project.)

_Hand-written — `gen_script_index.sh` leaves it alone._
