# byMUMmer — assembly-vs-assembly coverage / gap assessment (MUMmer)

Align one assembly to another with **MUMmer** (`nucmer`), keep the good collinear blocks,
then check how the query assembly covers the reference's **N gaps** and which genes sit in
covered vs uncovered gaps — i.e. how much a new assembly closes another's gaps.

Needs MUMmer (`nucmer`, `show-coords`, `delta-filter`, `mummerplot`) and the repo helpers
`join_coords.pl`, `deal_table.pl`, `deal_fasta.pl`. See `cmd_list` for the full example.

| script | role |
|--------|------|
| `stat1_mergeCov.pl` | merge the (manually curated) alignment blocks into per-region coverage |
| `stat2_chk_gapCover.pl` | tag each reference N-gap by whether the other assembly covers it (`ngs.nlis merged.stat`) |
| `stat3_addGeneTag.pl` | add gene tags to the gap-coverage table (`prot.gff3 coverTag`) |

---

## Workflow (see `cmd_list`)
```sh
nucmer --prefix=pb_ngs --minmatch 50 --mincluster 500 --breaklen 10000 ref.fa qry.fa
show-coords --rcl pb_ngs.delta > pb_ngs.coords
delta-filter -1 -i 95.0 -l 10000 -u 50 -o 10.0 pb_ngs.delta > pb_ngs.filter   # + mummerplot to eyeball

perl join_coords.pl pb_ngs.coords > pb_ngs.coords.jn        # join collinear blocks, then length/overlap filter
#   ... select good blocks, then MANUALLY curate -> manual_alignments.txt
perl stat1_mergeCov.pl manual_alignments.txt > manual_alignments.txt.stat

deal_fasta.pl -listSite '[nN]+' ref.scf.fa > ngs.nlis       # list the reference N gaps
perl stat2_chk_gapCover.pl ngs.nlis manual_alignments.txt.stat > ngs.nlis.coverTag
perl stat3_addGeneTag.pl ref.annot_prot.gff3 ngs.nlis.coverTag > ngs.nlis.coverTag.geneTag
```

_Hand-written — `gen_script_index.sh` leaves it alone._
