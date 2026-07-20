#!/bin/bash
# ============================================================================
# run-all.sh — reproduce EVERY claim about Alpoge's Jacobian-Conjecture
# counterexample in two independent computer algebra systems (Singular and
# PARI/GP), plus negative controls. No sympy anywhere. Absolute binary paths
# only, per house rules.
#
#   Usage:  bash /home/alex/code/jc/verify/cas/run-all.sh
# ============================================================================
set -u
SING=/usr/bin/Singular
GP=/usr/bin/gp
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 2

pass=0
fail=0

run() {
  # $1 = tool binary, $2 = script file, $3 = human label
  local tool="$1" script="$2" label="$3"
  echo "============================================================"
  echo ">>> $label"
  echo ">>> $tool $script"
  echo "------------------------------------------------------------"
  local t0 t1 out
  t0=$(date +%s.%N)
  if [ "$tool" = "$SING" ]; then
    out="$("$tool" -q "$script" 2>&1)"
  else
    out="$("$tool" -q "$script" 2>&1)"
  fi
  t1=$(date +%s.%N)
  echo "$out"
  # tally PASS/FAIL tokens
  local p f
  p=$(printf '%s\n' "$out" | grep -c ': PASS')
  f=$(printf '%s\n' "$out" | grep -c ': FAIL')
  pass=$((pass + p))
  fail=$((fail + f))
  awk -v a="$t0" -v b="$t1" -v p="$p" -v f="$f" \
    'BEGIN{printf(">>> wall: %.3fs   PASS:%d  FAIL:%d\n", b-a, p, f)}'
  echo
}

echo "############################################################"
echo "# SINGULAR — positive claims"
echo "############################################################"
run "$SING" claim1_det.sing            "CLAIM 1  Jacobian det == -2"
run "$SING" claim2_fiber.sing          "CLAIM 2  fiber over (-1/4,0,0): witnesses + vdim=3"
run "$SING" claim3_generic_fiber.sing  "CLAIM 3  generic fiber count == 3"
run "$SING" claim4_eliminants.sing     "CLAIM 4  x- and y-eliminants (Groebner elimination)"
run "$SING" claim5_discriminant.sing   "CLAIM 5  discriminant square factor (resultant)"
run "$SING" claim6_special_fibers.sing "CLAIM 6  empty fibers (unit ideal) + single point (vdim=1)"
run "$SING" claim7_equivariance.sing   "CLAIM 7  G_m equivariance, weights (-2,-1,+1)"
run "$SING" claim8_punt.sing           "CLAIM 8  Punt covering identity"

echo "############################################################"
echo "# SINGULAR — negative controls (breaks MUST be detected)"
echo "############################################################"
run "$SING" neg_controls.sing          "NEG CONTROLS (Singular)"

echo "############################################################"
echo "# PARI/GP — positive claims (independent CAS)"
echo "############################################################"
run "$GP" pari_claims.gp                "CLAIMS 1-8 (PARI/GP)"

echo "############################################################"
echo "# PARI/GP — negative controls"
echo "############################################################"
run "$GP" pari_negatives.gp             "NEG CONTROLS (PARI/GP)"

echo "############################################################"
echo "# TOTALS"
echo "############################################################"
echo "PASS tokens: $pass"
echo "FAIL tokens: $fail"
if [ "$fail" -eq 0 ]; then
  echo "RESULT: ALL CHECKS PASS (every claim reproduced in Singular + PARI; every negative control detected its break)."
  exit 0
else
  echo "RESULT: $fail FAILURE(S) — inspect output above."
  exit 1
fi
