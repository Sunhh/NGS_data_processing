#!/bin/bash
# Regression harness for deal_table.pl (characterization tests).
# Usage: ./run_tests.sh [update]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/../../deal_table.pl"
PMLIB="$HERE/../../MyPM"
GOLD="$HERE/golden"; MODE="${1:-run}"; cd "$HERE"; mkdir -p "$GOLD"
CASES=$(cat <<'EOF'
column02	-column 0,2 fixtures/t.tab
kick1	-column 1 -kick_col fixtures/t.tab
reverse	-reverse fixtures/t.tab
skip1	-skip 1 fixtures/t.tab
maxcol	-max_col 2 fixtures/t.tab
mincol	-min_col 2 fixtures/t.tab
colstat	-col_stat 2 fixtures/t.tab
coluniq	-col_uniq 1 fixtures/t.tab
colreps	-col_reps 1 fixtures/t.tab
colrepcount	-col_repCount 1 fixtures/t.tab
colsort	-col_sort 2 -col_sort_rule -1 fixtures/t.tab
transpose	-transpose fixtures/t.tab
fillNull	-fillNull NA fixtures/ragged.tab
cbind	-cbind fixtures/t.tab fixtures/t2.tab
combine	-combine fixtures/t.tab fixtures/t2.tab
uniqcolline	-UniqColLine 1 fixtures/t.tab
bestuniq	-best_uniqCol 1 -select_col 2 -select_rule 1 fixtures/t.tab
colhead	-col_head fixtures/t.tab
ksrch	-kSrch_idx fixtures/idx.tab -kSrch_idxCol 0 -kSrch_srcCol 0 fixtures/t.tab
dR2dN	-dR2dN fixtures/crlf.tab
trimEnd	-trimEndReturn fixtures/crlf.tab
EOF
)
pass=0; fail=0
while IFS=$'\t' read -r name args; do
  [ -z "$name" ] && continue
  out=$(perl -I"$PMLIB" "$SCRIPT" $args 2>/dev/null)
  if [ "$MODE" = update ]; then printf '%s\n' "$out" > "$GOLD/$name.out"; echo "  updated $name"
  else
    if diff -q <(printf '%s\n' "$out") "$GOLD/$name.out" >/dev/null 2>&1; then pass=$((pass+1))
    else fail=$((fail+1)); echo "  FAIL: $name"; diff <(printf '%s\n' "$out") "$GOLD/$name.out" 2>&1 | head -6 | sed 's/^/      /'; fi
  fi
done <<< "$CASES"
[ "$MODE" = update ] && { echo "golden regenerated ($(ls "$GOLD"|wc -l) files)."; exit 0; }
echo "=== deal_table regression: $pass passed, $fail failed (of $((pass+fail))) ==="
[ "$fail" = 0 ]
