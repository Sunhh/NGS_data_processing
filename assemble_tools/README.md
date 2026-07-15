# `assemble_tools` — genome assembly, scaffolding & QC toolbox

A grab-bag of scripts used across genome assembly: aligning two assemblies to
**join scaffolds**, editing **AGP / scaffold structure**, detecting and
**breaking mis-assemblies** from read coverage, **removing redundant** contigs,
and QC (BUSCO, QUAST, k-mer, LAI). Technology-specific helpers live in
subdirectories.

## Subdirectories

| dir | purpose |
|-----|---------|
| [`hifi_hic/`](hifi_hic/) | HiFi / Hi-C assembler output -> FASTA / AGP (hifiasm, HiCanu, 3D-DNA) |
| [`pacBio_tools/`](pacBio_tools/) | PacBio PBJelly gap-filling log parsing |
| [`bionano/`](bionano/) | BioNano optical-map `.xmap` helpers |
| [`map_ctg_to_scf/`](map_ctg_to_scf/) | contig -> scaffold placement (RagTag AGP) |
| [`classify_tools/`](classify_tools/) | contamination / organelle screening vs NCBI nt |
| [`kmer/`](kmer/) | k-mer counting (Jellyfish), GenomeScope summaries |
| [`high_tandem_repeat/`](high_tandem_repeat/) | tandem-repeat k-mer profiling |
| [`busco/`](busco/) | BUSCO `full_table` summaries |
| [`LAI/`](LAI/) | LTR Assembly Index plotting |

---

## Top-level scripts by task

### Join scaffolds by whole-genome alignment (LAST + Mugsy)
Align two assemblies, build a synteny MAF, and chain the linkage into merged
scaffolds. `list_run_last2scaf.pl` is the **driver** that wires the rest together.

| script | role |
|--------|------|
| `list_run_last2scaf.pl` | **driver**: LAST-align a list of assembly pairs, Mugsy WGA, then link+join scaffolds (set `$MUGSY_SRC` to your Mugsy install) |
| `run_last_to_scaffold.sh` | shell recipe for the LAST alignment step |
| `run_mugsy_MP.pl` | run Mugsy WGA in parallel (`-mugsy_src` / `-mugsy_exe` override the defaults) |
| `cmd_batch_for_mugsy` | batch command template for Mugsy |
| `clip_scaf_end.pl` | trim scaffold ends (N-runs) before aligning |
| `add_tag_to_fsa.pl` | prefix Ref/Qry tags onto FASTA IDs |
| `format_maf_forMugsy.pl`, `maf2fasta.pl`, `get_paired_maf.pl` | reshape MAF for the pipeline |
| `good_link_fromMaf.pl` | extract confident scaffold linkages from the synteny MAF |
| `grp_maf_byLinkage.pl` | group scaffolds by linkage |
| `join_link.pl` | join scaffolds along the linkage graph -> super-scaffolds |

### AGP / scaffold-structure editing
| script | role |
|--------|------|
| `scf2LG_to_AGP.pl` | linkage-group scaffold order -> AGP |
| `link_scf2chr.pl` | place scaffolds onto chromosomes |
| `link_seq_by_agp.pl` | build scaffold/chromosome FASTA from an AGP |
| `order_scf_1.pl` | order scaffolds by a position table |
| `extract_ctg_from_scaf.pl` | split scaffolds back into contigs at N-gaps |
| `mk_bed_from_agp.pl`, `lift_bed_jcviPlot.pl` | AGP -> BED / lift BED through an AGP (JCVI plots) |
| `add_1bp_ctg_toAGP_jcviPlot.pl` | pad tiny contigs into an AGP for JCVI plotting |
| `cnvt_loc_fromAGP_toAGP_forLoci.pl` | lift locus coordinates between two AGPs |
| `resize_N_in_agp.pl` | rewrite gap (N-run) sizes in an AGP |
| `estimate_gap_size.pl` | estimate gap sizes from paired-end spanning reads |
| `fill_SingleNgap.pl` | fill single-N gaps where flanks join cleanly |
| `filter_dropGap.pl` | drop gap-only / trivial records |

### Mis-assembly detection & breaking (read coverage)
| script | role |
|--------|------|
| `brk_fas_by_FRBadWind_for15kb.pl` | break contigs at windows with bad forward/reverse read pairing |
| `mk_wind_for_INScnt.pl` | make windows for insert-count scanning |
| `draw_aln_from_bam.pl` | draw alignments over a scaffold region from BAM |
| `get_frag_cov.pl` | fragment (insert) coverage from `samtools view -h in.bam` |
| `add_dep_to_chopInf.pl` | attach per-base depth to a chop-info table |
| `ctgBn6_to_scfCov.pl` | contig->scaffold coverage from a blastn (bn6) + AGP |
| `get_rep_loc_fromPileup.pl`, `loc_repeat.pl` | flag high-depth (repeat) locations from pileup |
| `estimate_gap_size.pl` | (see AGP) also uses PE coverage |

### Remove redundant sequence
| script | role |
|--------|------|
| `rmRed_byIdentCov_InCtg.pl` | drop short contigs covered by longer ones (megablast, identity+coverage) |
| `rm_redundant_loci.pl` | drop redundant loci from a nucleotide FASTA (blastn self-align) |
| `deduplicate_ncpu.pl` | QUAST-based dedup within & among many FASTAs (parallel, pan-genome) |

### QC / misc
| script | role |
|--------|------|
| `cnvt_quast_unaligned_info_to_tbl.pl` | QUAST unaligned-contig report -> table |
| `calc_ident_from_sam2pairwise.pl` | percent identity from `sam2pairwise` output |
| `slct_pe.pl` | select paired-end read subsets |
| `mk_hapmap_from_SNPtbl.pl` | SNP table -> HapMap (for map-based scaffolding) |
| `cnvt_loci_97_v1_to_vRILs.pl` | project-specific locus-coordinate conversion |
| `hist_plot_ins.R` | insert-size histogram |
| `plot_kmer_prop_along_chr.r` | k-mer proportion along chromosomes |
| `funcs_for_compare_genome.r` | R helpers for genome comparison plots |

---

_Hand-written navigation index (not the `gen_script_index.sh` autogen). External
tools: LAST, Mugsy, blastn/megablast (BLAST+), samtools, QUAST, BUSCO, Jellyfish,
GenomeScope, RagTag, hifiasm/HiCanu/3D-DNA, BioNano, LTR_retriever._
