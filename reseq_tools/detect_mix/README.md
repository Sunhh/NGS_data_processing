# detect_mix — introgression / admixture from group major alleles

Classify each individual's ancestry along the genome using group **major-allele (mjAL)**
patterns, then count the classification per window to reveal introgressed / mixed segments
(a run where an individual's genotypes match a different group's major alleles).

The input is a **mjAL table** — the major allele of each reference group at each site
(`chr, pos, grp1, grp2, …`, N when undefined), computed upstream from per-group genotypes.

| script | role |
|--------|------|
| `class_mjAL.pl` | mjAL table → per-site class list (the group major-allele pattern, i.e. which groups are diagnostic at each site) |
| `class_snpTab_by_mjAL.pl` | SNP table + class list → assign each individual's genotype at each site to a group class (`-class_mjAL`, `-colN_start 3`) |
| `class_grpCnt_in_mjTab.pl` | per-window count of each individual's group classification (`-ind2grp_list`, `-wind_length/-wind_step`) → introgression track |

---

## Workflow
```sh
# 0. mjAL table: major allele per group per site (chr pos grp1 grp2 ...), made upstream.

# 1. classify sites by the group major-allele pattern
#    -onlyDiag: keep only diagnostic sites (>=2 distinct group alleles)
perl class_mjAL.pl -onlyDiag set02.mjAL.CLV_CLM_CA_CC > set02.mjAL.class_list

# 2. assign each individual's genotype at each site to a group class
perl class_snpTab_by_mjAL.pl set02.snp.tab -class_mjAL set02.mjAL.class_list -colN_start 3 > set02.mjAL_class

# 3. per-window per-individual group counts (the introgression / mixture track)
perl class_grpCnt_in_mjTab.pl set02.mjAL_class -ind2grp_list grp_diff_list \
  -wind_length 100000 -wind_step 100000 > set02.mjAL_class_byWind
```
A window where an individual is repeatedly classified into a group other than its own
flags a candidate introgressed / admixed segment.

## Method notes & limitations
This is a fast, reference-free **diagnostic-allele** heuristic — good for a first-pass view
of ancestry blocks, but statistically limited; interpret candidate segments with care.
- Uses only each group's **major allele**, ignoring its frequency: a 51%-vs-99% difference
  counts the same. Prefer `-onlyDiag` (drops non-diagnostic sites); a frequency-difference
  filter would need per-group allele frequencies (not in this input format).
- A group whose major allele is **`N`** (undefined) is now **left out** of a site's class by
  default (it no longer wildcard-matches every allele and get double-assigned). Use
  `-nAsWildcard` to restore the old permissive behaviour.
- No probabilistic model / significance / recombination: it cannot separate real gene flow
  from shared ancestral polymorphism (ILS) or genotyping error.
- For rigorous introgression testing use **Dsuite** (D / f_d / f-branch, controls for shared
  ancestry); for probabilistic local ancestry with recombination use **RFMix / AncestryHMM /
  Loter / ELAI / MOSAIC**.


_Hand-written — `gen_script_index.sh` leaves it alone._
