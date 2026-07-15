# SV_detection — structural variants from whole-genome alignment

Call SNPs and structural variants (InDels, duplications, relocations, …) between a reference
and query genome from a **whole-genome alignment**, convert to VCF, add read-coverage support,
and annotate affected genes. The main path is **AnchorWave → NucDiff**; minimap2→paftools and
Assemblytics are alternatives.

External tools: `anchorwave`, `minimap2`, MUMmer (`nucmer`), **NucDiff** (patched copy in
[`nucdiff_modification/`](nucdiff_modification/)), **Assemblytics**
([`assemblytics_scripts/`](assemblytics_scripts/)), `last` (`maf-convert`), `picard`, `samtools`.
The `.sh` drivers and `cmd_list_example_to_detect_sv` are project run-scripts — edit the paths
at the top for your data.

## Pipeline (AnchorWave → NucDiff), see `callsv_AW_NDF.sh` / `cmd_list_example_to_detect_sv`
1. **Align**: `anchorwave proali` (CDS-anchored) → whole-genome collinear MAF; fix it with
   `../../software_fix/anchorwave/fix_awMAF.pl`.
2. **MAF → delta**: `maf-convert sam` → `sam2delta.py` (MAF/SAM helpers below prepare it).
3. **Call**: `nucdiff` on the delta → `*_ref_snps.gff` + `*_ref_struct.gff`.
4. **GFF → VCF**: `cnvt_ndfGff2vcf_snps.pl` / `cnvt_ndfGff2vcf_struct.pl`.
5. **Clean / support / annotate**: coverage, filtering, gene overlap (tables below).

## Scripts by role
**Alignment → delta (MAF / SAM / PAF / anchors)**
| script | does |
|--------|------|
| `sam2delta.py` | SAM → nucmer `.delta` (NucDiff input) |
| `fix_sam_cigarID.pl`, `restore_sam_position.pl`, `if_needIDfix.pl`, `join_samAln.pl` | fix/join alignment SAM records |
| `join_maf_blk.pl`, `rm_0span_maf.pl`, `filter_maf_by_tab.pl`, `get_shrt_or_ident_mafTab.pl`, `mmp2Aln_anchorsInMafTab.pl` | join / filter MAF blocks |
| `fmt_paf.pl`, `make_fakeCDS_fromPAF.pl` | minimap2 PAF → table / fake-CDS |
| `cnvt_anchors_to_tbl.pl`, `view_anchors.pl`, `find_nonOvlCDS.pl` | AnchorWave anchors |
| `align_2seq_by_query_segments.pl`, `get_qryLoc_by_refLoc_inBam.py` | pairwise realignment / coord lookup |

**NucDiff GFF → VCF**
| script | does |
|--------|------|
| `cnvt_ndfGff2vcf_snps.pl` | NucDiff SNP GFF → VCF |
| `cnvt_ndfGff2vcf_struct.pl` | NucDiff structural GFF → VCF (handles tandem-dup / relocation / N-gap edges; 2023 update) |
| `get_qry_ref_shared_var_nucdiff.pl` | variants supported by both ref- and qry-based NucDiff runs |
| `gff2bed.pl` | NucDiff GFF → BED (`-for_nucdiff`) |

**VCF cleanup / representation**
| script | does |
|--------|------|
| `cnvt_vcf_sub2insdel.pl`, `cnvt_vcf_sub2ins.pl` | recode long substitutions as INS/DEL |
| `remove_gap_var.pl`, `rmdup_fromNormVcf.pl`, `rmNvar_inVCF.pl` | drop gap-edge / duplicate / N variants |

**Read-coverage support & filtering** — `add_rdCov2vcf.pl` (add Q2R/R2Q depth), `filter_rdCov_vcf.pl`

**Select / summarise / annotate genes** — `get_sv_inVCF.pl`, `select_var.pl`, `summary_svs.pl`,
`get_sv_affected_genes.pl`, `add_geneAHRD.pl`

**Drivers** — `run_ndf.sh` (NucDiff), `run_mm2nucdiff.sh` / `run_mm2paftool.sh` (minimap2),
`callsv_AW_NDF.sh` (AnchorWave+NucDiff), `cmd_list_example_to_detect_sv`.

Bundled third-party (patched): `nucdiff_modification/` (NucDiff `*_revised`), `assemblytics_scripts/`.

_Hand-written — `gen_script_index.sh` leaves it alone._
