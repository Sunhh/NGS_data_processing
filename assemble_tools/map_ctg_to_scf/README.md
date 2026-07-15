# `assemble_tools/map_ctg_to_scf` — contig→scaffold placement (RagTag)

| script | role |
|--------|------|
| `stat_ragtag_agp.pl` | summarise a RagTag `ragtag.scaffold.agp` (per chromosome, by ID regex) |
| `simple_fill.pl` | place query contigs into reference scaffolds from a blastn (bn6) mapping |

```sh
perl stat_ragtag_agp.pl '^Cla97Chr' BC19/.../ragtag.scaffold.agp > stat.tsv
```
_Hand-written._
