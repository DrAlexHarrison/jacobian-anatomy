#!/bin/bash
# One-command reproduction of every computational claim in NOTE.md.
# Requires: python3 with sympy.  Exits nonzero on any failed certificate.
set -e
cd "$(dirname "$0")"
for f in verify/01_basic.py verify/02_anatomy.py verify/03_new_facts.py \
         verify/04_parametrization.py verify/05_presentations_rigidity.py \
         verify/06_identities.py verify/07_second_counterexample.py; do
  echo "=== $f"
  python3 "$f"
done
if command -v Singular >/dev/null 2>&1; then
  echo "=== verify/g4_fibers.sing (Singular)"
  OUT=$(Singular -q verify/g4_fibers.sing 2>&1)
  echo "$OUT"
  echo "$OUT" | grep -q "collision fiber vdim: 8" || { echo "FAIL g4 collision"; exit 1; }
  echo "$OUT" | grep -q "generic fiber vdim: 8"   || { echo "FAIL g4 generic"; exit 1; }
  echo "$OUT" | grep -q "5x8-104x4+95x2+4"        || { echo "FAIL g4 eliminant"; exit 1; }
else
  echo "(Singular not found — skipping verify/g4_fibers.sing)"
fi
echo
echo "ALL SUITES PASS."
echo "Lean certificate: run  'lake env lean JC.lean'  inside a mathlib4 checkout."
