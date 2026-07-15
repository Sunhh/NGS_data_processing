# NGS_data_processing — refactor log & optimization playbook

A running record of the MyPM library overhaul + dependent-script migration
(July 2026), plus a playbook for continuing to optimize the rest of this repo.
Start here next time.

---

## Working rules (read first)

- **Edit MyPM in place, here** (`MyPM/`). The old external working copy
  `/data/sunhh/temp/claude_code/MyPM` has been deleted; this repo is the single
  source of truth.
- **Never run networked git in this repo** (no `git push` / `pull` / `fetch` /
  `clone`). Local `git status` / `diff` / `log` for inspection only; commits are
  made manually.
- **MyPM has a regression test net** at `MyPM/t/`. After any MyPM change:
  `cd MyPM && prove -I. t/` (46 tests). Add assertions before refactoring a sub.
- **Backward compatibility is the constraint**: many `.pl` scripts depend on MyPM.
  Keep every public sub name/signature working; verify with the test net + a
  behavioral old-vs-new probe.
- **mathSunhh now requires its `mathSunhh/` subdirectory** (the loader `require`s
  the topical parts). Any copy/deploy of MyPM MUST include `mathSunhh/`.

---

## Status & resume here  (updated 2026-07-11)

> **DECISION (2026-07-11):** the owner has decided to leave `reseq_tools/snpTbl_stats.pl` as-is — do NOT spend further effort on it (the BioPerl-free `snpTbl_stats_v2.pl` remains available for anyone who wants it). It is the only MyPM-using script that fails to compile, and that is now intentional/accepted.

### Done so far (each script verified; test net in `tests/<name>/` unless noted otherwise)
| target | verify | key changes |
|---|---|---|
| **MyPM library** | `MyPM/t/` (46) | mathSunhh split into `mathSunhh/*.pm` + pure-functional (no OO); cross-module dedup; bareword calls package-qualified; PopGenSunhh (replaces Bio::PopGen); `parseCol` made public |
| **deal_fasta.pl** | `tests/deal_fasta/` (24) | GBK→UTF-8; dead code; `rcSeq`+`get_fasta_seq`→MyPM; fixed `-aa2cds` bug; parseCol |
| **deal_table.pl** | `tests/deal_table/` (21) | `tmsg`/`stop`→LogInforSunhh; parseCol |
| **temp/deal_gff3.pl** | `tests/deal_gff3/` (6) | `rev_comp`→fastaSunhh::rcSeq |
| **temp/get_cds_from_gff3.pl** | reused deal_gff3 fixtures | `rev_comp`→fastaSunhh::rcSeq (file is CRLF) |
| **deal_fastq.pl** | baseline-compare | double-mojibake encoding fix; `rcSeq`→fastaSunhh |
| **annot_tools/predictByAug_rna2genome.pl** | reader-equiv + `-c` (pipeline, not runnable) | `get_fasta_seq`→ReadInSeqSunhh; dead `use` removed |
| **assemble_tools/kmer/get_kmer_by_seq.pl** | `-w -c` clean | dead+shadowing `tsmsg` deleted |
| **assemble_tools/fill_SingleNgap.pl** | baseline-compare | `tsmsg`→LogInforSunhh |
| **solQ2phredQ.pl, ColLink.pl** | baseline-compare | dead-code only; kept local tmsg/parseCol (standalone — see policy) |
| **3 cmd scripts** | `-c` | dead `use lib` removed (run_cmd_in_batch.pl, cmd_ctrl/wait_for_pid.pl, cmd_ctrl/kill_by_ppid.pl) |
| **parseCol (cross-cutting)** | — | promoted to public `mathSunhh::parseCol`; ALL repo callers migrated; `_parseCol` kept as alias |
| **temp/detect_syn_dots.pl** | formula probe + `-c` | broken `Statistics::Regression` → inline closed-form OLS (verified pts (1,2)(2,4)(3,5)→slope1.5/intcpt0.667) |
| **evolution_tools/draw_syn_dotplot.pl** | renderer unit-test + `-c` | broken `DBD::Chart::Plot` → self-contained `SynDotPlot` SVG class (new/setOptions/setPoints/plot); out_fmt png→svg |

NOTE (2026-07-11): the repo's `.git/` directory has been DELETED by the owner — there is no git
history and no git-based rollback anymore. The working tree is the single source of truth; do not
run any git command here. (Earlier notes below about "review & commit" / "rollback via git" are moot.)

### Candidate sweep (local reimplementations of MyPM functions): COMPLETE for MyPM-users
Cleared categories — every script that ALREADY depends on MyPM: **get_fasta_seq, rev_comp,
rcSeq, tsmsg/tmsg**. Remaining local copies live only in STANDALONE scripts (no MyPM deps):
ColLink.pl, solQ2phredQ.pl, combine_DEGs.pl, extract_fq_by_list.pl,
assemble_tools/get_frag_cov.pl, assemble_tools/good_link_fromMaf.pl — left as-is by policy.
Re-scan anytime:
```
grep -rlnE '^sub (tmsg|tsmsg|stopErr|parseCol|rcSeq|rev_comp|get_fasta_seq)\b' --include='*.pl' .
```
**Policy:** delegate only when the script ALREADY uses MyPM AND behavior/format matches;
never force `use mathSunhh`/`fastaSunhh` (whole MyPM stack; mathSunhh needs CPAN
Statistics::Descriptive) onto a self-contained script.

### Next up — pick one (playbook)
> Dependency-hygiene (broken CPAN) status: Bio::PopGen ✓(v2), Statistics::Regression ✓, DBD::Chart::Plot ✓,
> Statistics::R ✓, Statistics::Descriptive::Weighted ✓(dead import), Statistics::Distributions ✓(chisqrprob).
> **Substantially DONE (2026-07-11).** NOTE: the earlier "COMPLETE" was premature — it only checked 4 named
> deps. A full load-test of EVERY imported module (see revised Diagnostic below) found 3 more broken-dep
> scripts; 2 are now fixed. Remaining: `snpTbl_stats.pl` (Bio::PopGen; v2 exists) and `converter_iprV4.pl`
> (legacy EBI iprscan-4.0 framework — NOT self-containable, left as-is).
1. **Dependency hygiene** — CPAN-replaceable broken deps are fixed (see status line above). Two loose
   ends remain, both by design: `reseq_tools/snpTbl_stats.pl` still imports Bio::PopGen (switch callers to
   the BioPerl-free `snpTbl_stats_v2.pl`, or port the fix in); and `annot_tools/iprscan/converter_iprV4.pl`
   needs the old EBI iprscan-4.0 custom framework (`Index::InterPro`, `Dispatcher::Tool::InterProScan`) —
   not self-containable, left as legacy.
2. **Hardcoded tool-paths** — ~76 `.pl` (mostly `run_*.pl`) hardcode machine exe paths; delicate
   (often intentional defaults) — needs per-file judgment.
3. **Per-subdir READMEs / index** for the ~590 `.pl`. **DONE 2026-07-11** (see "Per-directory README / master index" below; `gen_script_index.sh`).
4. **Grow test nets** for the most-used scripts.
Portability item — dead `use lib` sweep: **DONE**.

## What was done (July 2026)

### MyPM library
1. Removed dead machine-specific `use lib` lines (fileSunhh, LogInforSunhh).
2. Fixed plotSunhh's bogus `@EXPORT = qw(ins_calc)` (ins_calc lives in mathSunhh).
3. De-duplicated cross-module subs — `permutations`/`combinations` (SNP_tbl) and
   `_setHashFromArr` (gffSunhh) now delegate to the canonical copies in mathSunhh.
4. Wrote `MyPM/README.md` — module inventory, dependency layering, conventions.
5. Built the regression test net `MyPM/t/` (compile check + 45 unit assertions).
6. Split the 2574-line god-module `mathSunhh.pm` into a thin loader +
   `mathSunhh/*.pm` topical parts (Stats, Interval, Group, Windows, Encode, Color,
   IndexDB, …), all still `package mathSunhh` (one flat namespace).
7. Unified named-arg parsing: every module calls `&mathSunhh::_setHashFromArr(@_)`
   (removed the mathSunhh objects held only for this).
8. Decoupled internal bareword cross-calls (`&min` → `&mathSunhh::min`, etc.) so
   there is no hidden same-package coupling.
9. **Made mathSunhh pure-functional** (retired its OO side): removed
   `new`/`_initialize`, the `$ms_obj` singleton, and the dual-mode guards;
   `newNumber()` now uses a package-level pool (`$mathSunhh::_pkgObj`).
10. Added `MyPM/PopGenSunhh.pm` — self-contained population-genetics stats
    (pi, theta, segregating_sites_count, tajima_D) replacing `Bio::PopGen::Statistics`.

Each step was verified byte-equivalent to the previous state (test net + probes).

### Dependent scripts (this repo)
- **Migrated 34 `.pl` (118 call-sites)** from `$obj->method(...)` to
  `mathSunhh::method(...)` — required because mathSunhh is now function-only.
  Purely mechanical; `git diff` shows only those two patterns. These are tracked
  by git — **review and commit** when ready.
- Created **`reseq_tools/snpTbl_stats_v2.pl`** = BioPerl-free copy of
  `snpTbl_stats.pl` using PopGenSunhh. The original `snpTbl_stats.pl` is untouched
  (still imports the unsupported `Bio::PopGen::Statistics`).

### Deferred / open
- Validate `snpTbl_stats_v2.pl` against a historical `snpTbl_stats.pl` window
  (owner will do). Note: v2 excludes `N` from **all** stats; v1 (BioPerl + the old
  `_selfTheta`) counted N — expect differences only at N-containing sites.
- Commit the untracked files: `MyPM/PopGenSunhh.pm`, `reseq_tools/snpTbl_stats_v2.pl`.
- `snpTbl_stats.pl` still needs Bio::PopGen — switch callers to v2, or port the fix
  into it, when ready.
- The `*.pm_ori` snapshots and `SNP_tbl.pm_bk` were removed (2026-07-10);
  rollback for MyPM is via git history.

---

## Playbook: optimizing the rest of the repo (~590 `.pl`)

Apply the same discipline used on MyPM — measure the surface, preserve behavior,
verify — one topic/subdir at a time. Ideas, roughly high-value first:

1. **Portability sweep.** Grep all `.pl` for hardcoded absolute paths, dead
   `use lib`, and commented-out blocks; centralize discovery via `PERL5LIB`.
2. **Dependency hygiene.** Find scripts importing unmaintained CPAN modules
   (e.g. `Bio::PopGen`, old BioPerl). Replace with self-contained code —
   `PopGenSunhh.pm` is the template.
3. **De-duplicate across scripts.** Many `.pl` copy-paste the same helpers
   (option parsing, temp files, table I/O). Move recurring logic into MyPM.
4. **Consistency of MyPM usage.** Audit which MyPM subs each script uses
   (`grep -rn 'subname' --include='*.pl' .`); retire dead/legacy paths (e.g. the
   two `get_fasta_seq` readers).
5. **Per-subdir READMEs / index.** 590 scripts across topical dirs
   (assemble_tools, reseq_tools, annot_tools, evolution_tools, …) — a short index
   per subdir would help a lot.
6. **Grow the test net** to cover the most-used MyPM subs and a few key scripts.

### Handy commands
```bash
cd MyPM && prove -I. t/                         # run MyPM regression tests
perl -I MyPM -c some_script.pl                  # compile-check a script
grep -rn 'mathSunhh::' --include='*.pl' .       # who calls what
git -C . status                                 # review changes (never push/pull)
```

_Last updated: 2026-07-11._

---

## Per-script progress

### `deal_fasta.pl` (1936→1912 lines; 16-yr swiss-army fasta tool, ~40 modes)
- **2026-07-10 — Tier 1 done (zero behavioral risk, verified):**
  - Converted file GBK→UTF-8 (old comments on 14 lines were mojibake; all confirmed
    to be *inside comments* — code byte-identical). Now valid UTF-8, tool-friendly.
  - Removed dead code: `sub check_CDS` (called 0×, marked "never use again") and a
    commented-out `tsmsg` stub.
  - Verified: `perl -c` OK; 7 modes (baseCount, N50, attribute, frag, upper, cds2aa,
    listSite) byte-identical to a pre-change baseline.
- **Regression harness built (2026-07-10):** `tests/deal_fasta/` — `fixtures/`,
    `golden/`, and `run_tests.sh` (23 characterization cases covering the Tier-2
    dedup targets: rcSeq via frag_r/frag_c/listSite/cds2aa/drawByList; parseCol via
    replaceID/drawByList; plus chop_seq, cds2aa, and broad mode coverage). Run:
    `bash tests/deal_fasta/run_tests.sh` (or `... update` to re-baseline after an
    intentional change). Verified it catches regressions.
- **BUG FOUND (pre-existing, not from Tier 1):** `-aa2cds` mode is broken — it calls
    unqualified `&setup_codon_tbl()` (resolves to main::, undefined) instead of
    `&fastaSunhh::setup_codon_tbl()`; dies "Undefined subroutine" (line ~422; the
    call at line ~274 is commented out). One-word fix, but it is logic — do it in Tier 2
    (with a test case added). Mode excluded from the harness for now.
- **Tier 2 done/decided (2026-07-10):**
    - `rcSeq` (7 calls) -> delegated to `&fastaSunhh::rcSeq`, local copy removed (verified byte-identical across seqs/tags; test net 24/24).
    - Fixed `-aa2cds` bug (`&setup_codon_tbl()` -> `&fastaSunhh::setup_codon_tbl()`); re-added its test case.
    - `site_list`/`chop_seq`/`cds2aa`/`aa2cds`: mode-handlers that ALREADY call MyPM internally -> nothing to dedup.
    - `get_fasta_seq` (30 calls): DELEGATED 2026-07-10 to `ReadInSeqSunhh::get_fasta_seq` (user accepted fail-hard). Normal FASTA byte-identical (verified reader-vs-reader + test net 24/24); a header lacking an ID now aborts ("Lack of Key") instead of yielding an empty key.
    - `parseCol`: left (behavior differs from mathSunhh::_parseCol). `Disp_seq`: no MyPM equivalent, left.

- **Superseded Tier 2 notes:**
  - `rcSeq` (7 calls): tr-table byte-identical to `fastaSunhh::rcSeq` → safe to delegate.
  - `site_list`/`chop_seq`/`cds2aa`/`aa2cds`: overlap fastaSunhh; verify interfaces first.
  - `get_fasta_seq` (31 calls): ≈ `ReadInSeqSunhh::get_fasta_seq`; high value, high risk — defer, needs a regression harness.
  - `parseCol` (2 calls): **NOT** equivalent to `mathSunhh::_parseCol` (that one also
    handles negative columns / skips empties) — do not blindly swap.
  - Prerequisite: build a small regression harness (sample fasta through key modes) before touching logic.

### `deal_table.pl` (1380→1371 lines; ~40-mode table swiss-army tool)
- **2026-07-10 — done:**
  - Tier 1: (encoding) — CORRECTED 2026-07-11: header comment lines 5-8 were GBK mojibake
    (the earlier "clean" claim was wrong); converted GBK->UTF-8 (comments only, code byte-identical, test net 21/21) and no dead code
    (the grep-flagged `usage`/`col_sort` are false positives: `&usage;` no-parens,
    and `col_sort` is a `sort` comparator).
  - Built regression harness `tests/deal_table/` (21 cases: column/kick, reverse,
    skip, max/min, col_stat, uniq/reps/repCount, col_sort, transpose, fillNull,
    cbind, combine, UniqColLine, best_uniqCol, col_head, kSrch, dR2dN, trimEnd).
    Run: `bash tests/deal_table/run_tests.sh`. Verified it catches regressions.
  - Tier 2 dedup: local `tmsg`/`stop` were byte-identical reimplementations of
    `LogInforSunhh::tsmsg`/`stopErr` (already imported) → replaced 6 `&tmsg`+10
    `&stop` calls with `&tsmsg`/`&stopErr`, removed the two local subs. Verified:
    test net 21/21, and the stopErr error path works at runtime (timestamped + exit 1).
  - Left as-is: `parseCol` (differs from mathSunhh::_parseCol — no negative cols),
    `is_digital`, `loc_tbl2hash` (specific).

### Cross-script: `parseCol` consolidated into `mathSunhh::_parseCol` (2026-07-10)
- There were 3 near-identical `parseCol` implementations (deal_fasta, deal_table,
  mathSunhh::_parseCol). Verified a 3-way comparison: **identical on all normal
  column specs** (`0,2`, `1-3`, `5,3-1`, `0,2,5-7`, …). `mathSunhh::_parseCol` is a
  strict superset (also handles negative columns, skips empty segments, tolerates
  internal whitespace).
- Delegated both scripts to it: deal_fasta (2 calls) and deal_table (15 calls) now
  call `&mathSunhh::_parseCol(...)`; local `sub parseCol` removed from both;
  `use mathSunhh;` added to deal_table's top. `mathSunhh::_parseCol` itself unchanged.
- Verified: deal_table 21/21, deal_fasta 24/24, MyPM 46/46.
- Edge-case behavior changes (do NOT affect normal specs, so real usage unaffected):
  negative columns now accepted (were errors); empty/trailing-comma segments now
  skipped (deal_table previously ERRORED on a trailing comma — this is a robustness
  improvement; deal_fasta previously pushed ''); unparsable tags now use
  mathSunhh's timestamped stopErr.
- Note: `_parseCol` is underscore-private by convention but now used by scripts; could
  be promoted to a public `mathSunhh::parseCol` later if desired.

### `mathSunhh::_parseCol` promoted to public `mathSunhh::parseCol` (2026-07-10)
- Renamed `_parseCol` -> `parseCol` (in mathSunhh/Color.pm), kept
  `*_parseCol = \&parseCol` as a backward-compatible alias, added `parseCol` to the
  loader's `@EXPORT_OK`. Reason: `_parseCol` was already a de-facto public API used by
  ~8 repo scripts.
- Switched the two actively-refactored scripts (deal_fasta, deal_table) to the public
  `mathSunhh::parseCol`. The other ~6 scripts still calling `mathSunhh::_parseCol` keep
  working via the alias (verified they compile); they can be migrated later.
- **2026-07-10 follow-up:** migrated ALL remaining callers — every `.pl` in the repo now
  calls `mathSunhh::parseCol` (8 files updated, 0 residual `mathSunhh::_parseCol` calls).
  The `*_parseCol = \&parseCol` alias is kept in Color.pm as a safety net for any external
  callers. All changed scripts compile (snpTbl_stats.pl fails only on the unrelated
  pre-existing `use Bio::PopGen::Statistics`). deal_fasta 24/24, deal_table 21/21.
- Verified: parseCol==_parseCol (same sub), deal_table 21/21, deal_fasta 24/24, MyPM 46/46.

### `temp/deal_gff3.pl` (1802→1794 lines; gff3 dispatch tool, 32 subs)
- **2026-07-10 — done:**
  - Tier 1: clean encoding, no dead code (all `action_*` handlers dispatched).
    Script already integrates MyPM heavily (gffSunhh, fastaSunhh::rcSeq,
    mathSunhh::compare_number_list/switch_position/ovl_region/_setHashFromArr, fileSunhh loaders).
  - Built focused regression harness `tests/deal_gff3/` (6 cases): the key one is
    `seqret -extractFeat CDS` on a gff3 with BOTH a forward and a **reverse-strand**
    multi-CDS gene (exercises the rev_comp path), plus getLoc mRNA/CDS, listTopID,
    simpleSort, sort. Run: `bash tests/deal_gff3/run_tests.sh`. Catches regressions.
  - Tier 2 dedup: `rev_comp` (1 call, in action_seqret) reimplemented the exact tr
    table of `fastaSunhh::rcSeq` while the script *already used* rcSeq elsewhere →
    delegated the call to `map { my $t=$_; &fastaSunhh::rcSeq(\$t,'rc'); $t } @sub_seqs`,
    removed the local sub. Verified: rev_comp==map+rcSeq byte-identical, test net 6/6
    (reverse-strand CDS unchanged).
  - Left as-is: `_rmClose_blk`/`_rmOvlap_blk` (specialized gff block filters, not
    mergeLocBlk duplicates), and other gff-specific helpers.

### `annot_tools/predictByAug_rna2genome.pl` (499→468 lines; AUGUSTUS pipeline wrapper)
- **2026-07-10 — done:**
  - Delegated local `get_fasta_seq` (1 call, in splitMfasta) to `ReadInSeqSunhh::get_fasta_seq`
    (added `use ReadInSeqSunhh`, removed local sub). Verified reader-vs-reader byte-identical
    on a multi-contig genome fasta (splitMfasta only uses key+seq). Fail-hard accepted.
  - Removed dead `use IPC::Open3;` / `use Symbol;` (imported, never used; `exeCmd` comes from
    LogInforSunhh which loads its own IPC::Open3).
  - No runnable test net: the script drives augustus/tophat/samtools end-to-end (can't run
    standalone here). Verification = reader-equivalence + `perl -c`.

### `solQ2phredQ.pl` (50→26 lines; standalone Solexa→Phred qual converter)
- **2026-07-10 — done:** removed 24 lines of dead commented-out code (old `%h`-based
  impl + `sol2san` sub, superseded by the `tr///` one-liner). Output verified
  byte-identical to a pre-change baseline (small fastq).
- **Did NOT delegate `tmsg`:** it is a fully standalone script (only `use strict/warnings`);
  delegating to `LogInforSunhh::tsmsg` would force a MyPM dependency onto a self-contained
  utility AND change the stderr format (`$tt ` vs `[$tt]`). Net negative — left as-is.
  (Lesson for the candidate list: a local `tmsg`/`stop`/etc. is only worth delegating when
  the script already depends on MyPM and the behavior/format matches.)

### `ColLink.pl` (150→147 lines; standalone column-link/join tool)
- **2026-07-10 — reviewed:** Did NOT delegate its local `parseCol` — ColLink is
  self-contained (only `use strict; Getopt::Long`); delegating to `mathSunhh::parseCol`
  would drag in the whole MyPM stack incl. the CPAN dep Statistics::Descriptive, plus
  change error handling (die→stopErr). Kept local. Only removed 3 dead comment lines;
  output verified byte-identical (key-link test).

### `temp/get_cds_from_gff3.pl` (143→134 lines; extract spliced CDS from gff3)
- **2026-07-10 — done:** delegated local `rev_comp` (1 call) to `&fastaSunhh::rcSeq`
  (already `use fastaSunhh`), removed local sub. Same pattern as deal_gff3. Verified:
  compiles; CDS output (forward + reverse-strand gene) byte-identical to baseline
  after canonicalizing record order.
- Notes: this file has **CRLF** line endings (use `\r?\n` in regex edits); and its
  output record order is hash-nondeterministic (pre-existing, not touched).

### `deal_fastq.pl` (1037→1028 lines; fastq swiss-army tool)
- **2026-07-10 — done:**
  - Encoding fix: 2 comment lines were **double-encoded** (original GBK read as Latin-1
    then saved as UTF-8). Recovered with `iconv -f UTF-8 -t ISO-8859-1 | iconv -f GBK -t UTF-8`
    (verified: only those 2 lines change; recovers proper Chinese 全局变量 / 制作对应反向互补序列).
    → Reusable recipe for other "double-mojibake" files (differs from deal_fasta's plain GBK).
  - Tier 2: `rcSeq` (7 calls, byte-identical to `fastaSunhh::rcSeq` incl. in-place interface)
    → delegated (`use fastaSunhh` added, calls → `&fastaSunhh::rcSeq`, local removed).
    Verified: compiles; `-search both` and `-frag -frag_c` (rcSeq path) byte-identical to baseline.

### `assemble_tools/kmer/get_kmer_by_seq.pl` (133→129 lines)
- **2026-07-10 — done:** deleted the local `sub tsmsg` — it was **never called** AND
  shadowed the imported `LogInforSunhh::tsmsg` (byte-identical), triggering a
  `-w` "Subroutine tsmsg redefined" warning. Removal is behavior-neutral (dead) and
  clears the warning. `perl -w -c` now clean.
- Left `sub rc` (ATGC-only reverse-complement, 1 call): differs from `fastaSunhh::rcSeq`
  on IUPAC ambiguity codes (rcSeq complements N/R/Y…, `rc` only ATGC); delegating would
  also add a `use fastaSunhh` dep. Identical for pure-ATGC/N k-mers; left as-is.

### `assemble_tools/fill_SingleNgap.pl` (216→213 lines; single-N gap filler)
- **2026-07-10 — done:** local `sub tsmsg` (called 5×, byte-identical to
  `LogInforSunhh::tsmsg`) delegated: added `use LogInforSunhh;`, removed the local sub —
  the 5 calls now use the imported one. (Script uses SeqAlnSunhh, which already pulls in
  LogInforSunhh, so no new dep.) Verified: `-w -c` clean; STDOUT byte-identical to a
  baseline run, STDERR identical modulo timestamps.

### Candidate sweep status (2026-07-10)
All candidates that **already depend on MyPM** are now done:
- get_fasta_seq: deal_fasta.pl, predictByAug_rna2genome.pl — **category cleared**
- rev_comp: temp/deal_gff3.pl, temp/get_cds_from_gff3.pl — **category cleared**
- rcSeq: deal_fasta.pl, deal_fastq.pl — **category cleared**
- tsmsg/tmsg: get_kmer_by_seq.pl (dead, deleted), fill_SingleNgap.pl (delegated) — MyPM-users cleared
Remaining local reimplementations are all in **standalone scripts** (no MyPM deps) —
per the rule above, left as-is (delegating would add heavy/CPAN deps): ColLink.pl,
solQ2phredQ.pl, combine_DEGs.pl, extract_fq_by_list.pl, assemble_tools/get_frag_cov.pl,
assemble_tools/good_link_fromMaf.pl. Re-scan anytime with the grep in "Next candidates".

### Portability sweep — dead `use lib` (2026-07-10, playbook item 1)
- Removed dead `BEGIN { use lib "/abs/path" }` blocks (same machine-specific/missing
  paths as MyPM's, e.g. `/usr/local/share/perl5/`, the dead conda ForkManager path) from:
  run_cmd_in_batch.pl, cmd_ctrl/wait_for_pid.pl, cmd_ctrl/kill_by_ppid.pl. All find
  LogInforSunhh via PERL5LIB; all compile. Repo now has no hardcoded `use lib` outside MyPM.
- STILL OPEN (bigger, delicate): ~76 `.pl` hardcode machine tool-paths (mostly `run_*.pl`
  wrappers with default exe paths like /home/sunhh/... , /Data/Sunhh/...). These are often
  intentional defaults — do NOT blanket-change; needs per-file judgment. Deferred.

## Dependency hygiene (playbook item 1) — 2026-07-10

**Diagnostic** — which third-party deps actually load here:
- WORK (leave alone): Bio::SeqIO, Bio::SearchIO, Bio::DB::Fasta, Bio::Location::* (BioPerl core
  installed), Statistics::Descriptive, SVG, Parallel::ForkManager.
- BROKEN (❌ won't load → the script can't run): `Bio::PopGen::*`, `Statistics::R`,
  `DBD::Chart::Plot`, `Statistics::Regression`, `Statistics::Descriptive::Weighted`,
  `Statistics::Distributions`, `XML::Quote`, and the custom `Index::InterPro` /
  `Dispatcher::Tool::InterProScan` (old EBI iprscan-4.0 framework, not on this machine).
- Rule: only replace a dep that is actually broken/uninstallable; don't rip out working BioPerl.
- **2026-07-11 method upgrade:** the original list was hand-picked. Now load-test EVERY imported
  module and grep its importers — reproducible re-scan:
  `for m in $(grep -rhoE '^\s*use\s+[A-Z][A-Za-z0-9_:]+' --include=*.pl . | sed -E 's/^\s*use\s+//' | sort -u); do perl -I MyPM -M"$m" -e1 2>/dev/null || echo "BROKEN $m"; done`

**Broken-dep scripts & status:**
- `reseq_tools/snpTbl_stats.pl` (Bio::PopGen) — replacement exists: `snpTbl_stats_v2.pl` (PopGenSunhh). [done earlier]
- `assemble_tools/get_frag_cov.pl` (Statistics::Descriptive::Weighted) — **DONE 2026-07-11**: the
  `use Statistics::Descriptive::Weighted` was a DEAD import (script computes with its own local `sub mean`;
  no Weighted object ever created). But the missing module made the whole script un-loadable (BEGIN failed).
  Removed the one line → now compiles + runs end-to-end (verified on a tiny SAM). Standalone script (no MyPM),
  so nothing else touched.
- `evolution_tools/ortho_tools/run_positive_selection.pl` (Statistics::Distributions) — **DONE 2026-07-11**:
  single call `chisqrprob($freedom,$delta_LRT)` (LRT p-value). Added pure-Perl `mathSunhh::chisqrprob($df,$x)`
  = chi-square upper-tail p via regularized incomplete gamma (Numerical Recipes gammln/gser/gcf), in
  mathSunhh/Stats.pm + @EXPORT_OK; 5 unit tests (t/10_math.t → 54 total). `use Statistics::Distributions`
  → `use mathSunhh`; call → `&mathSunhh::chisqrprob`. Verified vs standard chi-square critical values and
  R pchisq(...,lower.tail=F) to 7 decimals. Full run drives PAML/codeml (not runnable standalone here).
- `annot_tools/iprscan/converter_iprV4.pl` (Index::InterPro, Dispatcher::Tool::InterProScan, XML::Quote) —
  **LEFT AS-IS (legacy, 2026-07-11)**: this is the ~2003 EBI InterProScan-4.0 output converter; it depends on
  that distribution's custom Perl framework (not CPAN, not in this repo). Not self-containable without vendoring
  the whole old framework. Flagged, not fixed.
- `temp/detect_syn_dots.pl` (Statistics::Regression) — **DONE 2026-07-10**: it used a simple
  2-param regression (`const`+`someX` → `theta`); replaced with inline closed-form OLS
  (y = fA + fB·x). Was `BEGIN failed` (un-loadable) → now compiles/loads. OLS == unweighted
  Statistics::Regression theta (verified formula: pts (1,2)(2,4)(3,5) → fB=1.5, fA=0.667).
- `enrich/scripts/enrichment_mine_fit.pl` (Statistics::R) — **DONE 2026-07-11**: R was used ONLY for
  `p.adjust(data[,7], method="BH")` (Benjamini-Hochberg FDR). Added pure-Perl `mathSunhh::p_adjust_BH(\@p)`
  (in mathSunhh/Stats.pm; added to @EXPORT_OK; 3 unit tests in t/10_math.t → 49 total) and rewrote the script to
  collect p-values in row order + write `_adjp` itself. `use Statistics::R` → `use mathSunhh`. Script was
  un-loadable (BEGIN failed) → now compiles + runs end-to-end with no R. Verified byte-exact vs R on the
  known example p.adjust(c(0.01,0.5,0.9,0.001,0.2)) = 0.025 0.625 0.9 0.005 0.333333, plus cap/monotonicity
  cases, plus a full end-to-end enrichment run (p-values pulled up correctly by the cummin step).
- `evolution_tools/draw_syn_dotplot.pl` (DBD::Chart::Plot) — **DONE 2026-07-10**: added a small
  self-contained `SynDotPlot` package (drop-in for the used API: new/setOptions/setPoints/plot) that
  emits SVG (no external dep); removed `use DBD::Chart::Plot`; default out_fmt png→svg. Script was
  un-loadable → now loads & compiles. Renderer unit-tested (line series→<polyline>, points→<circle>,
  nopoints→none, title/axis labels, colors). NOTE: output is now SVG (browser-viewable / convertible
  to PNG), not DBD's raster. Full end-to-end run needs MCScanX synteny input (not fabricated); verified
  via unit test + `perl -c`.


BioPerl-core scripts (still WORK, low priority): annot_tools/zff2augustus_gbk.pl (Bio::SeqIO/DB::Fasta
— writes GenBank, hard to replace), temp/cnvt_pairwise_to_tab.pl (Bio::SearchIO — parses blast pairwise).

## MyPM-user audit & further dedup (2026-07-11)

Focus: the scripts that depend on MyPM packages. Applied the same "don't trust a
hand-picked list — enumerate everything" lesson from dependency hygiene.

### Structural health check — CLEAN
- **394 `.pl` use ≥1 MyPM module.** Batch `perl -I MyPM -c` on all 394 → only
  `reseq_tools/snpTbl_stats.pl` fails (Bio::PopGen; `_v2` exists). The mathSunhh
  OO→functional refactor broke **no** script at compile time.
- **No residual mathSunhh OO usage**: 0 hits for `mathSunhh->new`, `$x = mathSunhh->…`,
  or removed subs (`new`/`_initialize`). Migration to `mathSunhh::foo(...)` was complete.
- **All 93 distinct `Module::sub` references resolve** to a defined sub (loaded every
  MyPM module, checked `$mod->can($sub)`). No dangling qualified calls.
- `@EXPORT` change check: plotSunhh no longer exports `ins_calc` — its only 2 users
  (`plot_syn.pl`, `plot_syn_bk.pl`) also `use mathSunhh`, and nobody calls bareword
  `ins_calc`. Safe.

### Local-reimplementation re-scan (supersedes the 7-name candidate grep)
Enumerated every local `sub` in the 394 scripts (873 defs) and intersected names with
MyPM's 225 sub names → 92 collisions. Triage (extract both bodies, comment/ws-normalized
diff): 67 are noise (`usage`, `new`). Non-noise verdicts:
- **IDENTICAL → delegated (DONE):** `aref_cols2tab` in `reseq_tools/cnvt_tools/cols2tab.pl`
  and `evolution_tools/vcf_tab/cols2vcfTab.pl` — byte-identical to `SNP_tbl::aref_cols2tab`,
  both already `use SNP_tbl`, body calls only stopErr/tsmsg/uc (no script-local helpers).
  Removed the local sub+POD (−74 lines each), qualified the 2 calls → `&SNP_tbl::aref_cols2tab`.
  Verified: compiles + **byte-identical** cols→tab output before/after (documented example).
- **Behavior-SUPERSET → delegated (DONE, user-approved 2026-07-11):** all 4 already use
  LogInforSunhh. `runCmd` (deduplicate_ncpu.pl, enrich/scripts/stat_goslim.pl,
  reseq_tools/gatk/revert_alnBam_to_uBam.pl) — removed local subs; runCmd is @EXPORTed so
  bareword `&runCmd` resolves to the canonical one (revert's was already behaviorally identical;
  the other two differed only in the failure-path stopErr message text, now unified to
  "[Err] Failed at CMD:"). `change_procN` (run_cmd_in_batch.pl) — not exported, so the 2 calls
  are qualified to `&LogInforSunhh::change_procN`; canonical adds a guard so an empty nprocF value
  keeps the previous max instead of setting it to 0 (otherwise identical). Verified: all `-w -c`
  clean; imported runCmd tested (success + failure exit 1); change_procN tested good/bad/empty.
- **Genuinely different (same name, different behavior) — LEFT AS-IS:** `geno2num` (×7,
  14–47-line diffs), `setup_windows` (×3), `tbl2seq` (×2, 177-line diff), `_readInAln` (×2),
  `permutations` (maskClose_in_1col.pl). Also `aa2cds`/`chop_seq`/`siteList` in
  deal_fasta/deal_fastq are mode-dispatchers that already call fastaSunhh internally (not dups).

## Per-directory README / master index (2026-07-11, playbook item 3) — DONE

Deliverable: a script index for the 587 `.pl` across 77 directories.
- **71 per-directory `README.md`** (dirs that lacked one): a table of every `.pl` with a
  one-line synopsis auto-extracted from the script's own `perl $0` / `Usage: $0` line
  (~11% expose none → "(no synopsis found)").
- **Top-level `SCRIPTS_INDEX.md`**: directory overview (all 77 dirs, counts, links) + inlined
  tables for the 6 dirs whose README is hand-written or top-level, so every script is indexed.
- **`gen_script_index.sh`**: reproducible, idempotent generator. Auto files carry an
  `<!-- AUTOGEN gen_script_index.sh -->` marker; the generator only (re)writes marked/missing
  files and never touches a README without the marker. Re-run any time: `bash gen_script_index.sh`.
- The 10 hand-written READMEs (README.md, reseq_tools/, temp/, enrich/, MyPM/, plotting/,
  reseq_tools/bsa/, evolution_tools/copy_number_var/, rnaseq_tools/map_to_{genome,transcriptome}/)
  are preserved verbatim.
- Lesson (again): a first cut keyed "hand-written" off git-tracked state, which broke after the
  auto files were committed; the marker approach is state-independent. To improve synopsis
  coverage, broaden the `Usage:`/`$0` regexes in the EXTRACT block.

## SAM-FLAG logic consolidation (2026-07-11)

Trigger: "can sam_filter.pl and sam_flag_chk.pl be merged?" Answer: they do different
jobs (a human-facing flag explainer vs a SAM-stream filter), so merging the CLIs is not
worthwhile — but BOTH (plus assemble_tools/get_frag_cov.pl) reimplement SAM-FLAG logic
that already lives canonically in **MyPM/SeqAlnSunhh.pm** (`sam_flag_infor()` = per-bit
decode; `mk_flag('keep'=>,'drop'=>)` = the keep/drop engine, already used by ~25 scripts).
So the real fix is delegation to SeqAlnSunhh, not merging.

- **MyPM/SeqAlnSunhh.pm**: added `sam_flag_table()` — accessor for the canonical
  `%infor_flag` (bit → [bitpos,hex,char,desc]), shared with sam_flag_infor(); no existing
  caller affected. Test net 54/54.
- **sam_flag_chk.pl** (MyPM-user already, via LogInforSunhh) → decodes via
  sam_flag_table()+sam_flag_infor(). Byte-identical for flags 0..2047; range extended to
  0..4095 so bit 11 (supplementary) is decoded instead of rejected.
- **sam_filter.pl** (was standalone) → keep/drop engine replaced by one mk_flag() call;
  %is_output rebuilt as {0..4095 => keep?1:0} to preserve the main path AND the -showNumber
  quirk exactly. 350→284 lines. Verified byte-identical over a synthetic all-4096-flag SAM
  across 12 rule sets + -showNumber + NM%/XT + no-filter error path.
- **assemble_tools/get_frag_cov.pl** (was standalone) → keep engine replaced by mk_flag(),
  restricted to 0..2047 to preserve original behaviour (supplementary still excluded).
  175→145 lines. Verified: kept-set identical for 0..2047 + byte-identical SAM output.

**Policy note:** sam_filter.pl and get_frag_cov.pl were previously standalone (zero MyPM
deps). Delegating adds the MyPM/SeqAlnSunhh stack — normally avoided by the standalone
policy, but **explicitly user-approved here** to get a single source of truth for SAM-FLAG
logic. `use SeqAlnSunhh` only imports olap_e2e_A2B, so it does not collide with those
scripts' local tsmsg/stopErr/mean.

## Runtime-resolution audit — undefined-sub calls (2026-07-11)

Beyond `perl -c` (which does NOT catch bareword calls to undefined subs), audited all 396
MyPM-using `.pl` for calls to subs that are neither locally defined, imported (@EXPORT or
explicit use-list), nor a Perl builtin. Tool: `audit_undef_subs.pl` (heuristic; strips
POD/comments/heredocs/strings/regexes, uses `prototype("CORE::$n")` for builtins and each
module's real @EXPORT).

**Result: ZERO real undefined-subroutine calls.** The scan raised 39 candidates; every one
was manually confirmed to be a regex/string fragment that looks like `word(` (e.g. `\t(\d+)`,
`m!^GT(:|$)!`, `Expect(?:\d+)?`), not an actual call. (`snpTbl_stats.pl` couldn't be audited
— it fails to load on Bio::PopGen — but it is the intentionally-abandoned script.)

Combined with the earlier structural audit (all compile except snpTbl_stats; no mathSunhh OO
residue; all 93 qualified Module::sub refs resolve), this gives strong confidence the MyPM
refactor introduced no latent breakage in its dependent scripts.

## Root-directory (top-level) optimization sweep (2026-07-11)

Tier-1 pass over the 21 top-level `.pl`. Fixed:
- **deal_table.pl** — GBK header comments (lines 5-8) → UTF-8 (log's earlier "encoding clean" claim was wrong). Comments only; test net 21/21.
- **deal_iprscan.pl** — removed dead `use Time::Piece` (only referenced in a commented-out line).
- **bp0_2_bp6.pl** — removed the dead commented-out local `sub openFH` (superseded by fileSunhh::openFH).

Verified clean (nothing to do): no dead MyPM `use` (checked incl. OO `$obj = Module->new` usage,
e.g. follow_mcscan's fastaSunhh IS used via `$fs_obj->save_seq_to_hash`); no other dead CPAN import;
no dead local subs; no broken deps; `deal_fastq.pl`'s `use lib "$pl_dir/MyPM"` is an intentional
self-locating path (not a dead machine path); only deal_table was non-UTF-8.

Deferred / left by policy: hardcoded default tool-paths (deal_augustus.pl etc. — playbook item 2,
per-file judgment); deal_bnx.pl's small commented-out alternative blocks (borderline documentation);
Tier-2 refactors such as bp0_2_bp6.pl's ~6×-duplicated "flush record" block (blocks differ subtly —
Lambda clears all %info vs others clear only send — and there is no regression net, so higher risk).

## annot_tools optimization sweep (2026-07-11)

> **DECISION (2026-07-11):** owner said to leave `iprscan/converter_iprV4.pl` alone (legacy EBI iprscan-4.0 framework) — do not spend further effort on it; it stays as the one non-loading annot_tools script, intentionally.

Tier-1 pass over the whole annot_tools tree (135 `.pl`). Fixed:
- **zff2augustus_gbk.pl** — removed dead `use File::Spec` (only reference was the use line;
  its BioPerl deps Bio::SeqIO / Bio::DB::Fasta / List::Util are all genuinely used).
- **intron2exex.pl** — removed dead `sub reversecomplement` (never called).
- **repAnno_tools/ch_gff_to_tab.pl** — removed dead `sub olapLen` and `sub minDist` (never called).

Verified clean: all 135 compile EXCEPT `iprscan/converter_iprV4.pl` (legacy EBI iprscan-4.0
framework — Index::InterPro / Dispatcher::Tool::InterProScan / XML::Quote — already known, left);
encoding all UTF-8; no dead MyPM `use`; no local reimplementations of MyPM helpers (candidate
sweep clean here); no other broken CPAN deps; undefined-sub audit = 0 real findings (its one flag,
`MAP` in intron2exex.pl, is a filehandle used as `print MAP (...)`, not a sub).

Left by policy: `augustus.accuracy_calculator.pl` defines a local `accuracy_calculator` that
duplicates `fromBraker::accuracy_calculator` (fromBraker's was likely extracted from it), but the
script is a 35-line standalone with no MyPM deps — delegating would force the MyPM stack onto it.

## assemble_tools optimization sweep (2026-07-11)

> **Follow-up bug fix (2026-07-11):** get_frag_cov.pl also had an inverted if/else in the per-scaffold flush loop that made every non-final scaffold output garbage (0,-1) instead of real coverage; fixed to match the final-scaffold loop (verified on a 2-scaffold SAM).

Tier-1 pass over the assemble_tools tree (64 `.pl`). Fixed:
- **classify_tools/cnt_In_bp.pl** — **real bug**: line 19 `print join("\t", $)."\n";` — `$)`
  parsed as the special GID var, leaving join()'s paren unclosed → syntax error; the script
  could not compile/run at all. Fixed to `print join("\t", $k1, $inCnt{$k1})."\n";` (output
  columns inferred from the accumulator logic; verified on the script's embedded example →
  ptg000012l 10825, ptg000099l 500). NOTE: adjust columns if a different format was intended.
- **ctgBn6_to_scfCov.pl** — removed dead `sub sep_ctgID`.
- **add_dep_to_chopInf.pl** — removed dead `sub idx_by_Pos` (use mathSunhh kept — still used).
- **run_mugsy_MP.pl** — removed dead `use Parallel::ForkManager` (never used).

Verified clean: all 64 compile; encoding all UTF-8; no broken CPAN deps; undefined-sub audit
= 0 (no flags at all); no clean dedup candidates (get_frag_cov.pl / good_link_fromMaf.pl keep
their local tsmsg/stopErr — standalone-style, format differs from LogInforSunhh's).

## rnaseq_tools optimization sweep (2026-07-11)

Tier-1 pass over the rnaseq_tools tree (43 `.pl`; subdirs graft/map_to_genome/coexp/
fromOthers/map_to_transcriptome all included). Fixed:
- **DEGtool_withSizeFactor.pl** — **real bug**: line 351 `keys $$total{$gene}` (keys on a
  reference) was experimental in Perl 5.14-5.22 and fatal since 5.24 → the script did not
  compile on this Perl 5.30 ("Experimental keys on scalar is now forbidden"). Fixed to
  `keys %{$$total{$gene}}` (same behavior). Now compiles.

Verified clean: all 43 compile (after the fix); encoding all UTF-8; no broken CPAN deps;
no dead subs; no dead imports (IO::File in DEGtool and SVG in draw_SNP are heavily used —
heuristic false alarms); undefined-sub audit = 0 real (its 2 flags in fromOthers/
run_TMM_scale_matrix.pl are false positives: `GetOptions` is exported despite the
`use Getopt::Long qw(:config ...)` config-list syntax, and `lib` came from print strings).

Left by policy: map_to_genome/combine_DEGs.pl keeps its local `parseCol` (standalone script,
only strict/warnings — delegating to mathSunhh::parseCol would drag in the whole MyPM stack).

## reseq_tools optimization sweep (2026-07-11)

Tier-1 pass over the reseq_tools tree (125 `.pl`; all 16 subdirs incl. cnvt_tools/gatk/xpclr/
slct_sweep/fst/... included). Fixed:
- **SNP_effect.pl, SNP_effect_edit.pl** — GBK Chinese comments → UTF-8 (were invalid UTF-8);
  code lines byte-identical (verified with grep -a), both compile.
- **fst/run_hierfstat.pl** — removed dead `sub run_fst_R_old` (superseded by run_fst_R).
- **xpclr/prepare_xpclr_input_wiGmP.pl** — removed dead `sub load_gmP` (script reads gmP from
  an input column instead).
- **cols2vcf.pl + cnvt_tools/cols2vcf.pl** — removed dead `sub geno2num` (superseded by the
  per-site %baseNum numbering in the main loop).
- **phylo_tools/replace_ID_in_nwk.pl** — removed dead `use Text::Balanced qw(extract_bracketed)`
  (never called).

Verified clean: all compile EXCEPT the abandoned snpTbl_stats.pl (Bio::PopGen); no other broken
deps; undefined-sub audit = 0; ForkManager/Statistics::Descriptive/SVG imports all genuinely used.

Left: `sort` comparators my_sort / my_sort_01 (not dead — heuristic false alarm); local geno2num
in cols2ped/cols2fstat/cnvt_tbl2fstat (USED and behaviorally different from SNP_tbl::geno2num),
setup_windows (differs from mathSunhh's), permutations, tbl2seq — genuinely different, kept;
extract_sites_by_list.pl's local openFH (standalone script — delegating would add fileSunhh dep).

## evolution_tools optimization sweep (2026-07-11)

Tier-1 pass over the evolution_tools tree (93 `.pl`; subdirs SV_detection/structure/ortho_tools/
vcf_tab/expansion_tools/compare_assemblies/... all included). Fixed:
- ortho_tools/run_positive_selection.pl — removed dead `sub print_std`.
- SV_detection/cnvt_ndfGff2vcf_struct.pl — removed dead `sub filter_array`.
- expansion_tools/03.replace_geneID_in_orthomcl.pl — removed dead `sub grpID`.
- expansion_tools/01.clean_nwk.pl — removed dead `use Text::Balanced qw(extract_bracketed)`.

Verified clean: all 93 compile (the earlier dep-hygiene fixes hold — detect_syn_dots.pl,
draw_syn_dotplot.pl, ortho_tools/run_positive_selection.pl); encoding all UTF-8; no broken deps;
undefined-sub audit = 0 real (its one flag, `GT` in compare_assemblies/cnvt_ploidy_inVCF.pl, is
the `m!^GT(:|$)!` regex, not a sub); SVG / ForkManager / File::Copy / Statistics::Descriptive
imports all genuinely used.

Left: draw_syn_dotplot.pl's local `_readInAln` (defined AND used at line 72; the script does not
use mcsSunhh, so the name-collision with mcsSunhh::_readInAln is coincidental — not a dedup).

## Remaining directories — one-pass sweep (2026-07-11) — all CLEAN

Swept the remaining 18 dirs (108 `.pl`) in one pass: cmd_ctrl, enrich/scripts,
file_type_based/{Proc_Reads,Proc_Sam}, GBS, gm_tools, log_tools, pcr_tools,
project/watermelon_pan_phaseI, relate_loc_1to1, sample_scripts, self_interest,
site_search, software_fix/anchorwave, temp/{,forRonan,scripts,temp_fix_gff3}.

Result: **nothing to fix.** All 108 compile; encoding all UTF-8; no broken CPAN deps; no
dead subs; no dead imports (7 candidates — File::Which, FileHandle, IO::Handle, CPAN,
Bio::SearchIO, LWP::Simple — all verified genuinely used); no local reimplementations of MyPM
helpers; undefined-sub audit = 0 real (one flag, `t` in project/.../summary_covBp_byRdDep.pl,
is a `\t(\d+)` regex fragment). (enrich/scripts/enrichment_mine_fit.pl was already fixed earlier
in the dependency-hygiene work.)

### Whole-repo Tier-1 sweep: COMPLETE
Every directory has now had the Tier-1 pass (compile / encoding / broken-dep / dead-code /
dead-import / MyPM-dedup / undefined-sub). Directories swept: root, annot_tools, assemble_tools,
rnaseq_tools, reseq_tools, evolution_tools, + these 18. Intentionally-left non-loading scripts:
reseq_tools/snpTbl_stats.pl (Bio::PopGen) and annot_tools/iprscan/converter_iprV4.pl (legacy EBI).

## Tier-2: duplicate-script survey (2026-07-11)

9 basenames exist in two locations each, and every pair is **byte-identical** (md5 match):
- INTENTIONAL ARCHIVE (leave): project/watermelon_pan_phaseI/{cnt_TEIPRacc, list_IPRacc,
  pipe_for_functional_annotation, ret_maker_abinit_gff3}.pl — this dir is a frozen project
  snapshot (26 files incl. README.md + cmd_list); the copies are deliberate.
- Candidate consolidations (owner's call — destructive, and some may be deliberate convenience
  copies): reseq_tools/{cols2ped,cols2vcf}.pl == reseq_tools/cnvt_tools/{same} (cnvt_tools likely
  canonical); reseq_tools/cnvt_tools/cols2fstat.pl == reseq_tools/fst/cols2fstat.pl;
  gm_tools/{cnt_allele_withPop,cvt_snp_to_itayNeed}.pl == evolution_tools/vcf_tab/{same}.
No other script references these by path (only READMEs / SCRIPTS_INDEX.md).
Nothing deleted — awaiting owner decision.

> **CORRECTION (2026-07-11, owner):** the md5 survey above followed symlinks, so it
> conflated symlinks with real copies. Actual state (checked with `-type l` + inode):
> - **2 are SYMLINKS, not duplicates** (already single-source, the ideal arrangement):
>   project/watermelon_pan_phaseI/pipe_for_functional_annotation.pl -> annot_tools/ahrd/... ,
>   project/watermelon_pan_phaseI/ret_maker_abinit_gff3.pl -> annot_tools/maker/... .
> - The other **7 are real, separate byte-identical copies** (distinct inodes, not hardlinks):
>   cnt_allele_withPop.pl, cnt_TEIPRacc.pl, cols2fstat.pl, cols2ped.pl, cols2vcf.pl,
>   cvt_snp_to_itayNeed.pl, list_IPRacc.pl — left per owner.
> The whole repo has exactly 2 `.pl` symlinks (no directory symlinks). Where a real-copy
> pair was edited earlier (cols2vcf.pl geno2num removal) both sides were changed and remain
> byte-identical. NOTE for future: use `find -type l` / inode, not md5, to detect duplicates.

## Tier-2 follow-ups (2026-07-11)

- **get_frag_cov.pl output bug** — fixed (see note under the assemble_tools sweep): the
  per-scaffold flush loop's inverted if/else made every non-final scaffold print 0,-1.
- **deal_fastq.pl siteList redefine** — `use fastaSunhh` imported fastaSunhh::siteList (OO
  method) which deal_fastq's local `sub siteList ($$$)` then redefined ("Subroutine siteList
  redefined" under -w). Changed to `use fastaSunhh ()` (rcSeq is called qualified; the local
  siteList is the intended one). -search output byte-identical before/after.
- **Repo-wide `perl -w -c` scan** — the only remaining notable warnings are in the vendored
  third-party `annot_tools/repAnno_tools/ProtExcluder1.1/` scripts (main::usage / main::score
  "used only once") — left untouched as external code.

### Systematic bug-hunting: exhausted
Compile errors, encoding, broken deps, dead code/imports, MyPM-dedup, undefined-sub calls, and
-w redefine/used-once warnings have all been swept repo-wide. Real bugs found & fixed this
session: cnt_In_bp.pl (syntax), DEGtool_withSizeFactor.pl (keys-on-scalar), get_frag_cov.pl
(inverted output). Remaining work is judgment-heavy Tier-2 (e.g. bp0_2_bp6.pl's duplicated
flush block — needs a characterization test) or the delicate hardcoded-tool-path sweep.

## Portability: in-repo hardcoded tool paths -> FindBin $REPO (2026-07-11, playbook item 2)

~74 `.pl` hardcode absolute machine paths; ~50 point INTO this repo as
`/home/Sunhh/tools/github/NGS_data_processing/<rel>` (works here only via a
`/home/Sunhh -> sunhh` symlink; breaks on any other checkout). External-tool/data paths
(trimmomatic, bwa db, iprscan install, GO obo, GAF...) are genuine machine-specific
overridable defaults and are LEFT alone.

**Done (22 scripts):** replaced in-repo `//=` / `my $x=` defaults with a portable repo
root computed from the script's own location:
`use FindBin; (my $REPO = $FindBin::RealBin) =~ s{(/NGS_data_processing)(/.*)?$}{$1};`
then `... "$REPO/<rel>"`. The $REPO line goes right after the first `use` (before any use
of it). Verified: all compile AND every `$REPO/<rel>` default resolves to an existing file.
(generate_dataset.pl was the pilot; +21 in bccc4bf.)

**Left / flagged for a separate decision:**
- **6 stale defaults  14 FIXED (6a334fb):** these pointed at moved files (broken today), corrected to the real location AND portablized via $REPO  14 graft(3)->rnaseq_tools/map_to_genome/, structure(2)->evolution_tools/structure/, stat_goslim pl_obo2tab->enrich/scripts/ + ref_goslim->enrich/example_data/ (ref_bgobo is a user-provided external GO obo, left). All compile; targets verified to exist.
  resolve even today): enrich/scripts/stat_goslim.pl (enrich/... should be enrich/scripts?),
  evolution_tools/structure/{00.run_Structure,prepare_structure_input}.pl (structure/ is under
  evolution_tools/), rnaseq_tools/graft/{find_transmit_step1,find_transmit_stepX1,
  simple_pipe_find_graft_rd_SE}.pl (runHisat2_with2pass.pl / fix_NHnum.pl are under
  map_to_genome/), file_type_based/Proc_Reads/run_trimmomatic.pl (Proc_Reads/polyAT_adp.fa).
  Fixing these means correcting WHERE the default points — needs per-case confirmation.
- **~4 scripts already using `$HOME/tools/github/NGS_data_processing/...`** (rm_geneFrag,
  run_repClass, run_spaln_prot2genom, run_mugsy_MP) — already env-relative; left.
- **12 "tricky" scripts** — in-repo path appears inside a `runCmd("perl /path ...")` string
  (6) or a gatk config-heredoc table (6, e.g. pipe_gatk*.pl `pl_getSam  /home/Sunhh/...`);
  need careful per-form handling. Not done yet.

## Portability sweep — COMPLETE (2026-07-12)

All in-repo hardcoded tool paths handled. Final tally:
- **28 scripts** now resolve in-repo sibling-script/data defaults via a FindBin `$REPO`
  (`(my $REPO=$FindBin::RealBin)=~s{(/NGS_data_processing)(/.*)?$}{$1}`): 22 simple `//=`/assign
  (6db2afd,bccc4bf) + 6 command-string (8f536d2). run_trimmomatic's ILLUMINACLIP default too (e688dcd).
- **6 stale defaults fixed** (pointed at moved files — broken today): graft×3, structure×2,
  stat_goslim (6a334fb).
- **6 gatk example-config heredocs** (`<<'CFG'`, non-interpolating templates the user copies+edits):
  can't use $REPO; corrected the one stale path Proc_Sam/->file_type_based/Proc_Sam/ so the example
  is accurate; other example paths verified correct (e688dcd).
- **Left by design:** external tool/data defaults (trimmomatic install, bwa db, iprscan, GO obo/GAF,
  the ref_bgobo GO obo not shipped in-repo); scripts already using `$HOME`/`~` for the repo path
  (run_mugsy_MP, run_repClass, run_spaln_prot2genom, rm_geneFrag, list_run_last2scaf — already
  env-relative); and inline `# /path` doc-comments (cosmetic; the code defaults beside them are ported).

Verification for every change: `perl -c` compiles AND each `$REPO/<rel>` (or corrected absolute)
target verified to exist. No end-to-end pipeline runs (external tools) — compile + path-resolution only.

---

## 2026-07-15 — gm_tools ABH consolidation (owner decision on the tracked duplicates)
Resolved two items from the 2026-07-11 cross-dir duplicate list:
- **Retired** `gm_tools/{cnvt_snp_to_itayNeed, cvt_snp_to_itayNeed, cnvt_snp_to_itayNeed_cnt}.pl`.
  Folded their capability into `gm_tools/rnaseq_snp_tools.pl`: new `-out_itayCsv <pref>` on
  `-task get_abh` (per-chromosome R/qtl CSVs) and `-out_itayCnt <file>` on `-task cnt_pop_allele`
  (compact `chr pos Total A_N B_N AF%`). Smoke-tested: both new outputs are byte-identical to the
  retired scripts on clean sites. (get_abh already used the same `SNP_tbl` classifier as cvt.)
- **Dropped** `gm_tools/cnt_allele_withPop.pl`; the byte-identical `evolution_tools/vcf_tab/cnt_allele_withPop.pl`
  is kept as the vcf.tab-utility home.
- **Still present (owner kept):** `evolution_tools/vcf_tab/cvt_snp_to_itayNeed.pl` — byte-identical to the
  now-retired gm_tools copy and also reproducible via `rnaseq_snp_tools.pl -task get_abh -out_itayCsv`;
  left in vcf_tab per the vcf.tab-utility grouping (retire later if undesired).
