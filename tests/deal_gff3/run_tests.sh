#!/bin/bash
# Regression harness for temp/deal_gff3.pl (characterization tests).
# Usage: ./run_tests.sh [update]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/../../temp/deal_gff3.pl"
PMLIB="$HERE/../../MyPM"
GOLD="$HERE/golden"; MODE="${1:-run}"; cd "$HERE"; mkdir -p "$GOLD"
CASES=$(cat <<'EOF'
seqret_cds	-inGff fixtures/in.gff3 -scfFa fixtures/genome.fa -seqret -extractFeat CDS
getLoc_mRNA	-inGff fixtures/in.gff3 -getLoc mRNA
getLoc_CDS	-inGff fixtures/in.gff3 -getLoc CDS
listTopID	-inGff fixtures/in.gff3 -listTopID
simpleSort	-inGff fixtures/in.gff3 -simpleSort
sort	-inGff fixtures/in.gff3 -sort
EOF
)
pass=0; fail=0
while IFS=$'\t' read -r name args; do
  [ -z "$name" ] && continue
  out=$(perl -I"$PMLIB" "$SCRIPT" $args 2>/dev/null)
  if [ "$MODE" = update ]; then printf '%s\n' "$out" > "$GOLD/$name.out"; echo "  updated $name"
  else
    if diff -q <(printf '%s\n' "$out") "$GOLD/$name.out" >/dev/null 2>&1; then pass=$((pass+1))
    else fail=$((fail+1)); echo "  FAIL: $name"; diff <(printf '%s\n' "$out") "$GOLD/$name.out" 2>&1 | head -8 | sed 's/^/      /'; fi
  fi
done <<< "$CASES"
[ "$MODE" = update ] && { echo "golden regenerated ($(ls "$GOLD"|wc -l))."; exit 0; }
echo "=== deal_gff3 regression: $pass passed, $fail failed (of $((pass+fail))) ==="
[ "$fail" = 0 ]
