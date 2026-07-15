#!/bin/bash
# Regression harness for deal_fasta.pl (characterization tests).
# Golden outputs capture the CURRENT behavior; a diff means behavior changed.
# Usage:
#   ./run_tests.sh          run all cases, compare to golden/
#   ./run_tests.sh update   (re)generate golden/ from the current script
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/../../deal_fasta.pl"
PMLIB="$HERE/../../MyPM"
GOLD="$HERE/golden"
MODE="${1:-run}"
cd "$HERE"
mkdir -p "$GOLD"

# name <TAB> args (input files under fixtures/)
CASES=$(cat <<'EOF'
baseCount	-baseCount fixtures/t.fa
N50	-N50 fixtures/t.fa
attribute	-attribute key:len:GC:GCnum:AG fixtures/t.fa
frag	-frag 3-12 fixtures/t.fa
frag_r	-frag 3-12 -frag_r fixtures/t.fa
frag_c	-frag 3-12 -frag_c fixtures/t.fa
upper	-upper fixtures/t.fa
lower	-lower fixtures/t.fa
listSite	-listSite ATG -listBoth fixtures/t.fa
cds2aa	-cds2aa fixtures/cds.fa
aa2cds	-aa2cds fixtures/cds.fa fixtures/prot.fa
chop_seq	-chop_seq -chop_len 10 -chop_step 5 -chop_min 3 fixtures/t.fa
rmDefinition	-rmDefinition fixtures/t.fa
uniqSeq	-uniqSeq fixtures/dup.fa
uniqSeq_bySeq	-uniqSeq_bySeq fixtures/dup.fa
rna2dna	-rna2dna fixtures/rna.fa
keep_len	-keep_len 10-40 fixtures/t.fa
rmTailXN	-rmTailXN fixtures/tailN.fa
fq2fa	-fq2fa fixtures/t.fq
reorderByList	-reorderByList fixtures/order.txt fixtures/t.fa
replaceID	-replaceID -replaceIDlist fixtures/id_map.tab -replaceIDcol 0,1 fixtures/t.fa
drawByList	-drawByList -drawList fixtures/regions.tab -drawLcol 0,1,2,3,4 fixtures/t.fa
res	-res seq fixtures/t.fa
sample	-sample 1-2 fixtures/t.fa
EOF
)

pass=0; fail=0
while IFS=$'\t' read -r name args; do
  [ -z "$name" ] && continue
  out=$(perl -I"$PMLIB" "$SCRIPT" $args 2>/dev/null)
  if [ "$MODE" = update ]; then
    printf '%s\n' "$out" > "$GOLD/$name.out"; echo "  updated $name"
  else
    if diff -q <(printf '%s\n' "$out") "$GOLD/$name.out" >/dev/null 2>&1; then
      pass=$((pass+1))
    else
      fail=$((fail+1)); echo "  FAIL: $name"
      diff <(printf '%s\n' "$out") "$GOLD/$name.out" 2>&1 | head -6 | sed 's/^/      /'
    fi
  fi
done <<< "$CASES"

[ "$MODE" = update ] && { echo "golden regenerated ($(ls "$GOLD" | wc -l) files)."; exit 0; }
echo "=== deal_fasta regression: $pass passed, $fail failed (of $((pass+fail))) ==="
[ "$fail" = 0 ]
