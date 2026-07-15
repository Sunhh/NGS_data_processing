# compare_assemblies — compare genomes / assemblies

Assess one assembly against another and pull comparative features (overlaps, ploidy,
duplications) between genomes.

| item | does |
|------|------|
| [`byMUMmer/`](byMUMmer/) | MUMmer assembly-vs-assembly alignment → N-gap coverage + genes-in-gaps (how well one assembly closes another's gaps) |
| `cnt_ovl_fromBeds.pl` | count how many input BED files overlap each interval (`in_1.bed in_2.bed ... > ovl.cnt.bed`; uses bedtools) |
| `cnvt_ploidy_inVCF.pl` | recode VCF ploidy (e.g. 4x → 2x by randomly drawing 2 of 4 alleles) for cross-ploidy comparison |
| `mcscanTab_to_dupTxt.pl` | MCScanX self-self collinearity table → a duplicated-gene-pair list (WGD / segmental duplication) |

`mcscanTab_to_dupTxt.pl` consumes the MCScanX output produced via
`../ortho_tools/cnvt_gffJnLoc_to_mcscanx.pl`; synteny plotting is in the parent [`../`](../).

_Hand-written — `gen_script_index.sh` leaves it alone._
