#!/usr/bin/env bash
# Run run_dft.py over a list of stage-1 job dirs (one per line), low priority, one job at a time.
# Usage: run_queue.sh LIST OUTDIR [extra run_dft.py args]   (env: OMP_NUM_THREADS, PYSCF_MAX_MEMORY)
set -u
list=$1; out=$2; shift 2
here=$(cd "$(dirname "$0")" && pwd)
py=${PY:-$HOME/micromamba/envs/proton-relay/bin/python}
mkdir -p "$out"
while read -r d; do
  [ -z "$d" ] && continue
  echo "[$(date '+%F %T')] start $d"
  nice -n 19 ionice -c3 "$py" "$here/run_dft.py" "$d" --out "$out" "$@"
  rc=$?
  echo "[$(date '+%F %T')] end $d rc=$rc"
done < "$list"
echo "[$(date '+%F %T')] QUEUE DONE"
touch "$out/QUEUE_DONE.$(basename "$list")"
