#!/usr/bin/env bash
# Per-plot marginal ΔTmax sensitivity over the 400 cLHS plots, 12 chunks parallel.
# Reuses the Shapley 1111 baseline; runs 6 perturbations/plot. Resumable.
set -u
cd /home/corroyez/Documents/z_Example_rmusica_31012025
K=12
mkdir -p /tmp/perplot_logs
pids=()
for c in $(seq 1 $K); do
  Rscript c1_sensitivity_perplot_chunk.R "$c" "$K" > "/tmp/perplot_logs/chunk_${c}.log" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} chunks: ${pids[*]}"
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "ALL PERPLOT CHUNKS DONE (failures=$fail)"
ls out_files/Chapter3/nc_sensitivity_perplot | wc -l
cat out_files/Chapter3/tables/sensitivity_perplot/part_*.csv | grep -vc "^pid"
