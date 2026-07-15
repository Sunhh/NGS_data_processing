# `Proc_Reads` — FASTQ read processing

QC, cleaning, and mapping wrappers for raw Illumina reads. Thin Perl wrappers
build and run the external commands; a few `.sh` files are project recipes kept
as examples.

## Cleaning / QC
| script | role |
|--------|------|
| `run_fastqC.pl` | run FastQC on one or more FASTQs (`-inFq ... -outDir`) |
| `run_trimmomatic.pl` | **canonical Trimmomatic wrapper** — PE and SE in one batch tool (`-fq1_pref/-fq2_pref` prefix+suffix lists). Defaults are the historical trim combos; adapter fastas are in-package and parameter-specifiable (`-adp_PE`/`-adp_SE`/`-adp_polyAT`), jar defaults to the vendored `trimmomatic/trimmomatic.jar`, and `-bgzip` compresses outputs |
| `run_rmRRNA.pl` (`run_rmRRNA.sh`) | drop rRNA reads (SortMeRNA) over an input list |
| `srch_barcode.pl` | tally inline barcodes from `in_R1.fastq` |
| `run_ndupB.sh` | de-duplicate PCR duplicates (project recipe) |
| `chk_INS_byChlo.sh` | estimate insert size by mapping to a chloroplast reference (recipe) |

Adapter/contaminant FASTAs used by the trimming steps live here: `adaptors.fa`,
`illumina_adapters.fa`, `polyAT_adp.fa`. The vendored `trimmomatic/` is a
custom-built Trimmomatic (modified `SlidingWindowTrimmer`) plus its jar.

Poly-A/T tail removal is done by the R helper `using_subfunc.R` (a function
library); `using_subfunc.R.rm_polyAT.R` and `...rm_polyAT_useRight.R` are the two
trimming variants, `using_subfunc.R.cmd{,.R}` are runnable command examples.

## Mapping wrappers
| script | aligner |
|--------|---------|
| `run_bwaAln.pl` (`run_bwaAln.sh`) | `bwa aln`/`sampe` PE pipeline with step control (`-step2do aln_r1/aln_r2/sampe/bam_sort/...`) |
| `run_bowtie.pl` / `run_bowtie2.pl` | Bowtie / Bowtie2 (`-para_bwt` passthrough) |
| `run_tophat2.pl` | TopHat2 spliced RNA-seq mapping |

_Hand-written. Retired the superseded single-pair Trimmomatic wrappers
`run_trimmoPE.pl` / `run_trimmoSE.pl` and their project recipes
`cleanPE_byTrimmo.sh` / `cleanSE_byTrimmo.sh` (all hardcoded to non-portable
cluster paths; covered by `run_trimmomatic.pl`), plus the `using_subfunc.R.bk`
backup._
