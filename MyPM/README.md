# MyPM — Perl5 module library (Honghe Sun)

A personal bioinformatics Perl5 toolkit. Modules are named `<topic>Sunhh.pm`
and share a small foundation (`LogInforSunhh`, `mathSunhh`, `fileSunhh`).
All 14 modules compile under perl >= 5.30 and form a clean dependency DAG
(no cycles).

---

## 1. Installation / discovery

These are plain `.pm` files (not a CPAN dist). A script finds them by having
their directory on Perl's `@INC`. The recommended mechanism is `PERL5LIB`:

```bash
# in ~/.bashrc (use an idempotent, dedup-guarded append; see note below)
export PERL5LIB="/path/to/MyPM:$PERL5LIB"
```

Canonical live location on the author's machine:
`/home/sunhh/tools/github/NGS_data_processing/MyPM`

> Do **not** hardcode `use lib "/abs/path"` inside these modules (it bakes one
> machine's layout into the library). If a script must self-locate the library,
> put `use FindBin; use lib "$FindBin::RealBin/../MyPM";` in the *script*.
>
> If `.bashrc` is re-sourced per shell (screen/tmux windows), guard the append
> so directories are not duplicated in `PERL5LIB`.

### External (CPAN) dependencies

Only three are non-core; install them into your Perl environment (e.g. `~/perl5`):

| CPAN module            | Required by      | Used for |
|------------------------|------------------|----------|
| `Parallel::ForkManager`| LogInforSunhh    | `get_pm`, `change_procN` (multi-process) |
| `File::Which`          | fileSunhh        | `_which` |
| `File::Copy::Recursive`| fileSunhh        | `_dircopy` |
| `Statistics::Descriptive` | mathSunhh     | `ins_calc` and other stats |

Everything else (`IPC::Open3`, `Symbol`, `File::Copy`, `Cwd`, `File::Spec`,
`File::Path`, `File::Basename`, `Scalar::Util`) ships with core Perl.

---

## 2. Dependency layering

```
LogInforSunhh                 (foundation: logging, command exec, fork mgmt)
   └─ mathSunhh               (stats, interval/location math, grouping, encoding, color)
        └─ fileSunhh          (I/O, temp files, compression, table loaders)
             ├─ fastaSunhh    (FASTA / codon / translation)
             ├─ gffSunhh      (GFF3 read/write, overlap)
             ├─ SeqAlnSunhh   (SAM/BAM, BWA/TopHat, CIGAR, RBH/cscore)
             ├─ SNP_tbl       (SNP table manipulation & format conversion)
             └─ mcsSunhh      (MCScanX collinearity input parsing)
plotSunhh, ReadInSeqSunhh, ReadInAlnSunhh, ConfigSunhh, fromBraker, wm97Sunhh
   (leaf/utility modules; depend only on the foundation)
```

---

## 3. Module reference

Legend — Paradigm: **F** = functional (call `pkg::sub()` or import),
**OO** = object (`My->new`), **H** = hybrid.

| Module | Paradigm | Purpose | Default `@EXPORT` |
|--------|:--------:|---------|-------------------|
| **LogInforSunhh** | F | Timestamped logging, command execution, Parallel::ForkManager helpers | `tsmsg stopErr exeCmd exeCmd_1cmd runCmd` |
| **mathSunhh** | F | Statistics, interval/location arithmetic, grouping/graph, base-N encoding, RGB color, location-index DB | `ins_calc` |
| **fileSunhh** | F | File handles, temp file/dir, (de)compression, table/AGP/blast6 loaders, file splitting | `openFH renameByPat isSkipLine splitL wantLine wantLineC` |
| **fastaSunhh** | OO | FASTA reading, reverse-complement, codon tables, translation, 4dTv | `siteList` |
| **gffSunhh** | OO | GFF3 parse/read/write, feature typing, overlap between GFF sets | *(none)* |
| **SeqAlnSunhh** | OO | Read mapping (BWA/TopHat2), SAM/BAM parsing, CIGAR ops, SAM flags, e2e overlap, RBH/cscore | `olap_e2e_A2B` |
| **SNP_tbl** | OO | SNP-table filtering, genotype encode/decode, ts/tv, format export (fasta/nex/meg/ped/structure/hap/csv) | *(none)* |
| **mcsSunhh** | F | Parse MCScanX `.gff`/`.aln`/table, map chromosome position to plot coordinate | *(none)* |
| **plotSunhh** | F | HTML/hex/RGB color-name conversion; heat-map color (delegates to mathSunhh) | *(none)* |
| **ReadInSeqSunhh** | F | Streaming FASTA reader (positional-arg, legacy) | `get_fasta_seq` |
| **ReadInAlnSunhh** | F | MAF multiple-alignment reader / coordinate normalization | `readMAF splitMafSline` (+`normMAFloc` OK) |
| **ConfigSunhh** | OO | Read/write `key value` config files with `__var__` interpolation | *(none)* |
| **fromBraker** | F | AUGUSTUS accuracy parsing & config editing (from BRAKER) | *(none)* |
| **wm97Sunhh** | OO | Watermelon WM97 chromosome-ID ↔ number naming | *(none)* |
| **PopGenSunhh** | F | Population-genetics stats: nucleotide diversity (pi), Watterson theta, segregating sites, Tajima's D (self-contained; replaces Bio::PopGen::Statistics) | *(none)* |

---

## 4. Public subroutines by module

Subs prefixed with `_` are internal/private helpers (listed only where useful).

### LogInforSunhh
`tsmsg` · `stopErr` · `exeCmd` · `exeCmd_1cmd` · `runCmd` · `run` · `timed_out` ·
`usage` · `change_procN` · `get_pm` · `get_pid`

### mathSunhh (grouped)

> mathSunhh.pm is a thin loader; the subs live in `mathSunhh/*.pm` parts (all
> declare `package mathSunhh;`, so the namespace is one flat 'mathSunhh'). The themed groups
> below map to files: Stats->Stats.pm, Array/combinatorics->ArrayComb.pm, Object
> number->ObjNum.pm, Windows->Windows.pm, Interval->Interval.pm, Grouping->Group.pm,
> index DB->IndexDB.pm, Encoding->Encode.pm, Color->Color.pm, Param helper->Param.pm
> (constructor new/_initialize stays in mathSunhh.pm). **The `mathSunhh/` directory
> must travel with `mathSunhh.pm`.**
- **Stats:** `min` `max` `minmax` `log10` `ins_calc` · (private `_mean` `_sum` `_median`)
- **Array / combinatorics:** `repArr` `dvd_array` `permutations` `combinations`
  `randSlct_num` `create_randNum`
- **Object number mgmt:** `new` `newNumber` `offspringArray`
- **Windows:** `setup_windows` `map_windows`
- **Interval / location:** `ovl_region` `compare_number_list` `mergeLocBlk`
  `sep_loc2_by_loc1_multiLoc2` `sep_loc2_by_loc1_singleLoc2` `map_loc_byReference`
  `insert_SEloc_toSEArray` `switch_position` `transfer_position`
- **Grouping / graph:** `complete_pair_for_grouping` `divide_group`
  `divide_group_fromHash` `divide_group_fromArray` `extract_group`
  `extract_group_fromHash` `extract_group_fromArray`
- **Location/number index DB:** `index_SEloc` `index_numbers` `ret_sortedRealLoc`
  `ret_sortedNum` (+ a large set of `_`-prefixed internal indexing helpers)
- **Encoding:** `_Encode_deciN` `_Decode_deciN` `_decimal_to_hexa` `_hexa_to_decimal`
- **Color / columns:** `cnvt_to_rgb` `get_xy_byScale` `parseCol` (column-spec parser like `0,2,5-7`; public, in @EXPORT_OK; `_parseCol` kept as a legacy alias)
- **Param helper:** `_setHashFromArr` `_addHash`

### fileSunhh
`openFH` · `renameByPat` · `write2file` · `load_bn6File` · `load_tabFile` ·
`load_agpFile` · `reverse_agpHash` · `splitL` · `wantLine` · `wantLineC` ·
`isSkipLine` · `log_section` · `new_tmp_dir` · `new_tmp_file` · `dvd_file` ·
`iCompressFile` · `oCompressFile` · (thin File::* wrappers: `_which` `_copy`
`_move` `_dircopy` `_rmtree` `_abs_path` `_basename` `_dirname` `_catfile`)

### fastaSunhh
`new` · `save_seq_to_hash` · `get_fasta_seq` · `rcSeq` · `siteList` · `chop_seq` ·
`setup_codon_tbl` · `bbb2aa` · `aa2cds` · `aa2cds_1seq` · `get_4d_codon` · `cnt_4dtv`

### gffSunhh
`new` · `read_gff3File` · `write_gff3File` · `parse_line` · `ovl_between_gff3` ·
(internal: `_featType2Num` `_getAttrHash` `_allChild` `_allTopParent`
`_addParentID` `_setTypeHash` `_setHashFromArr` …)

### SeqAlnSunhh
Mapping: `aln_tophat2` `bwaAln` `bwaPE` `bwaSE` `chk_index` `version_samtools` ·
SAM/BAM: `openSam` `not_uniqBest` `sam_line2hash` `sam_hash_addKey`
`print_sam_lines` `cnt_sam_mismatch` · SAM flags: `mk_flag` `sam_flag_infor` ·
CIGAR: `cigar_str2array` `cigar_array2str` `cigar_array2len` `parseCigar`
`trim_pos` `trim_cigar_arr` `trim_cigar_str_bothEnd` · Homology:
`olap_e2e_A2B` `rbh_byBp6` `cscore`

### SNP_tbl
Object: `new` `newSubObj` `readTbl` `writeTbl` · Filter: `rm_lmiss` `rm_noVar`
`rm_multiVar` `max2Allele` · Genotype: `SingleChar` `SingleCharData` `get_tv`
`is_tstv` `cnt_maf` `cnt_genotype` `geno2num` `dna_d2b` `dna_b2d`
`get_diploid_d2b` `A2AA` `AA2array` · Cross/allele: `tab_allele`
`tab_allele_to_genotype` `tab_class_PP_al` `tab_class_off_al` `aref_cols2tab` ·
Export: `tbl2seq` `tbl2fas` `tbl2nex` `tbl2meg` `tbl2ped` `tbl2structure`
`tbl2hap` `tbl2illumina` `file_tbl2csv` · Misc: `guessChrNum` `chk_VarInArray`
`rmClose_idx`

### mcsSunhh
`chrP_to_plotP` · (internal readers: `_readInMcsGff` `_readInAln` `_readInAlnTbl`
`_readInChrLis`)

### plotSunhh
`_rgb_color` · `cnvt_to_rgb`

### ReadInSeqSunhh
`get_fasta_seq`  *(positional args `($fh, $has_head)`; legacy — see notes)*

### ReadInAlnSunhh
`readMAF` · `splitMafSline` · `normMAFloc`

### ConfigSunhh
`new` · `getConfig` · `writeConfig`

### fromBraker
`accuracy_calculator` · `setParInConfig`

### wm97Sunhh
`new` · `number_to_chrID` · `chrID_to_number`

### PopGenSunhh
`pi` · `theta` · `segregating_sites_count` · `tajima_D`  *(data model: `@inds` = list of `{marker=>[alleles]}`; 'N' excluded from per-site sample size. Replaces Bio::PopGen::Statistics; used by `reseq_tools/snpTbl_stats_v2.pl`.)*

---

## 5. Conventions & notes

- **Canonical ownership (de-duplicated).** These algorithms live in **mathSunhh**;
  same-named subs elsewhere are thin wrappers that delegate to it:
  - `permutations`, `combinations` — `SNP_tbl::*` delegate to `mathSunhh::*`.
  - `_setHashFromArr` — the single copy lives in mathSunhh; **all** modules parse
    their named args by calling the plain function `&mathSunhh::_setHashFromArr(@_)`
    (no module keeps a mathSunhh object just for this, and gffSunhh no longer has
    its own copy).
  - `cnvt_to_rgb` — `plotSunhh::cnvt_to_rgb` delegates to `mathSunhh::cnvt_to_rgb`.
- **Two FASTA readers (intentionally kept separate).**
  - `fastaSunhh::get_fasta_seq('faFh'=>.., 'faFile'=>.., 'has_head'=>..)` —
    OO, named args. **Preferred** for new code.
  - `ReadInSeqSunhh::get_fasta_seq($fh, $has_head)` — positional, exported.
    **Legacy**; kept for existing scripts.
- **Export policy.** Existing `@EXPORT` lists are a de-facto public API and are
  kept as-is. **New** functions should go to `@EXPORT_OK` (opt-in) rather than
  `@EXPORT`, to avoid polluting caller namespaces.
- **mathSunhh calling convention.** mathSunhh is **function-only**: every sub is
  called as `mathSunhh::foo(...)` (internal cross-calls are package-qualified too).
  It has no constructor and no objects. `newNumber()` keeps state in a package-level
  pool `$mathSunhh::_pkgObj` (one sequence per process).
- **Private subs.** A leading `_` marks an internal helper; treat it as unstable.
- The pre-refactor `*.pm_ori` snapshots and `SNP_tbl.pm_bk` were removed (2026-07-10)
  once verified; rollback for MyPM is via git history.

---

## 6. Change log & remaining ideas

For the full dated history of this refactor **and ideas for optimizing the rest of
the `NGS_data_processing/` repo**, see `../REFACTOR_LOG.md`.

**Still open / deferred:**
- Expand the `t/` regression suite to cover more subs.
- Validate `snpTbl_stats_v2.pl` (PopGenSunhh) against a historical `snpTbl_stats.pl` window (deferred; note v2 excludes N from all 4 stats, v1 did not).
- Commit the new/untracked files (`MyPM/PopGenSunhh.pm`, `reseq_tools/snpTbl_stats_v2.pl`) after review.

**Completed:**

- Expand the `t/` regression suite (started — see section 7) to cover more subs before larger refactors.
- ~~Split the large mathSunhh into cohesive units~~ **DONE** — mathSunhh.pm is now
  a loader over `mathSunhh/*.pm` parts (single `package mathSunhh` namespace, so all
  `mathSunhh::foo()` and `$obj->foo()` calls are unchanged). Verified byte-equivalent
  to the pre-split monolith across 25+ functions and the full test net.
- ~~Unify the two OO "param object" patterns~~ **DONE** — every module now parses
  named args via `&mathSunhh::_setHashFromArr(@_)`. Removed: the mathSunhh objects
  held *solely* for param parsing (fastaSunhh $mathObj, SeqAlnSunhh $ms, SNP_tbl
  $ms_obj) and gffSunhh's own _setHashFromArr wrapper. (In **Phase 2** below, gffSunhh's mathSunhh object was also removed: newNumber/
  offspringArray are now called as `mathSunhh::...` functions, newNumber using a
  package-level number pool.)
- ~~Retire mathSunhh's OO side (make it pure-functional)~~ **DONE (Phase 2)** —
  removed constructor/dual-mode/singleton; migrated 34 external `.pl` (118 call-sites)
  from `$obj->method` to `mathSunhh::method`; newNumber -> package-level pool.

---

## 7. Running the test net

A starter regression suite lives in `t/` (uses core `Test::More`; finds the
modules via `FindBin`). Run all of it with:

```bash
prove -I. t/          # from the MyPM/ directory
```

Current coverage:

| File | Covers |
|------|--------|
| `t/00_compile.t` | all 14 modules load (`require_ok`) |
| `t/10_math.t`    | mathSunhh: min/max/minmax/log10, permutations/combinations, _setHashFromArr, ovl_region, mergeLocBlk, repArr |
| `t/20_file.t`    | fileSunhh: splitL, isSkipLine, openFH roundtrip, load_tabFile |
| `t/30_fasta.t`   | fastaSunhh: rcSeq, siteList, bbb2aa, get_fasta_seq |
| `t/40_snp.t`     | SNP_tbl: dna_d2b/dna_b2d, combinations delegation |
| `t/50_readseq.t` | ReadInSeqSunhh: get_fasta_seq (legacy positional reader) |

These are **characterization tests**: expected values were captured from the
current implementation, so they lock in today's behavior. Re-run `prove` after
any change; a red test means behavior shifted. Add assertions here before
refactoring a sub, then confirm green afterward.

> Gotcha locked in by tests: `fileSunhh::splitL()` chomps its 2nd argument in
> place, so it must be given a modifiable variable (not a literal).
