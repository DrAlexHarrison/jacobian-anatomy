#!/bin/bash
# P2-EQUIVARIANT production sweep: 3x Groebner unit-certificate per system.
cd "$(dirname "$0")/.."
JOBS=${1:-jobs}; RES=${2:-results}; TMO=${3:-1750}
mkdir -p "$RES"
ls $JOBS/f*.sing | sort | xargs -P 8 -I{} bash -c '
  base=$(basename {} .sing)
  timeout '"$TMO"' /usr/bin/Singular -q {} > '"$RES"'/$base.out 2>&1
  echo "$base exit=$?" >> '"$RES"'/sweep.log
'
echo "SWEEP DONE $(date)" >> "$RES"/sweep.log
