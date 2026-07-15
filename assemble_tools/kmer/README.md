# `assemble_tools/kmer` — k-mer / genome-size utilities

k-mer counting queries (Jellyfish) and genome-size / heterozygosity summaries.

| script | role |
|--------|------|
| `get_kmer_by_seq.pl` | look up k-mer counts for a sequence in a Jellyfish `-jf_db` |
| `get_seq_by_kmer.pl` | grow a sequence from a seed k-mer walking a Jellyfish `-jf_db` |
| `get_kmer_by_seq_summary.pl` | summarise per-sequence k-mer count output |
| `cnvt_qualFa_to_wind_avgTbl.pl` | windowed average of a `.qual` fasta (window_size window_step) |
| `extract_genomescope_summary.pl` | pull genome size / heterozygosity from GenomeScope `summary.txt` |

_Hand-written._
