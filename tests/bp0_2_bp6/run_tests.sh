#!/bin/bash
# Characterization tests for bp0_2_bp6.pl (BLAST bn0 -> bn6 / -snp / -geno).
# Exercises all flush paths (Score=/> /Query=/Lambda) + listSNP. Run: bash run_tests.sh
# Re-baseline after an intentional change: bash run_tests.sh update
cd "$(dirname "$0")"
ROOT=../..; PL=$ROOT/bp0_2_bp6.pl; FX=fixtures; GD=golden
run () { # mode-args...
  case "$1" in
    bn6)  perl -I "$ROOT/MyPM" "$PL" -in "$FX/test.bn0" -out /dev/stdout 2>/dev/null ;;
    snp)  perl -I "$ROOT/MyPM" "$PL" -in "$FX/test.bn0" -snp  -bn6 "$FX/restrict.bn6" 2>/dev/null ;;
    geno) perl -I "$ROOT/MyPM" "$PL" -in "$FX/test.bn0" -geno -bn6 "$FX/restrict.bn6" 2>/dev/null ;;
  esac
}
if [ "$1" = update ]; then for m in bn6 snp geno; do run $m > "$GD/$m.txt"; done; echo "golden updated"; exit 0; fi
pass=0; fail=0
for m in bn6 snp geno; do
  if diff -q <(run $m) "$GD/$m.txt" >/dev/null; then pass=$((pass+1)); else echo "FAIL: $m"; diff <(run $m) "$GD/$m.txt" | head; fail=$((fail+1)); fi
done
echo "=== bp0_2_bp6: $pass passed, $fail failed (of $((pass+fail))) ==="
[ "$fail" -eq 0 ]
