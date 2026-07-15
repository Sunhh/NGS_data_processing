# `assemble_tools/busco` — BUSCO result summaries

Summarise a **BUSCO** run's `full_table.tsv` and clean up its scratch files.

| script | role |
|--------|------|
| `summary_busco_full_table.pl` | count Complete/Duplicated/Fragmented/Missing from `full_table` |
| `geneCopyN_busco_full_table.pl` | per-BUSCO copy number from `full_table.tsv` |
| `rm_busco_intermediate_files.pl` | delete a BUSCO output dir's bulky intermediate files |

```sh
perl summary_busco_full_table.pl full_table.tsv > full_table.cnt
perl geneCopyN_busco_full_table.pl full_table.tsv > full_table.tsv.cnt
```
_Hand-written._
