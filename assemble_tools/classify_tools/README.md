# `assemble_tools/classify_tools` — contamination / organelle screening

Screen assembled contigs against **NCBI nt** (blastn) and classify each region as
**In** (target taxon, keep) or **Ex** (contaminant/organelle/rDNA, drop), then
extract the regions to remove. Long contigs are segmented before blast so low
`-max_target_seqs`/`-max_hsps` don't lose sensitivity, and hit coordinates are
mapped back to the original contigs.

| script | role | step |
|--------|------|------|
| `temporary_prepare_rmCont.pl` | split over-long contigs by `max_length` before blast | 1 |
| `run_seg_blastn.pl` | segment query, blastn each segment vs `db`, lift positions back to origin coords | 2 |
| `line_bn6_query.pl` | collapse a bn6 table to one line per query | 3 |
| `classify_region_byBn6.pl` | classify each aligned region In/Ex from taxonomy-annotated bn6 | 4 |
| `get_Ex_region.pl` | emit separate/joined Ex (external) region lists | 5 |
| `recog_organelle_rDNA_from_classJn.pl` | tag chloroplast/mito/rDNA/bacteria classes | 5 |
| `cnt_In_bp.pl` | count In (kept) bp per contig | report |

```sh
perl run_seg_blastn.pl out_prefix 10000 to_rm_contamination.fa nt_db
perl line_bn6_query.pl out_prefix.bn6 > out_prefix.bn6.1query1line
perl classify_region_byBn6.pl ...            # -> In/Ex per region
perl get_Ex_region.pl P1R02scaf.nt_bn6.MiORGN_join 1>sep_ex_lis 2>joined_ex_lis
```
_Hand-written. blastn `-outfmt 6` must include taxonomy columns (staxids/sscinames/sskingdoms)._
