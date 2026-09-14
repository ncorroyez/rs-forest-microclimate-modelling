#!/usr/bin/env bash
# Regen global-baseline Shapley (legacy binary) — 12 chunks parallel, wait all.
# Full regen: the Blois global nc were deleted, so MuSICA runs all 6400.
set -u
cd /home/corroyez/Documents/z_Example_rmusica_31012025
K=12
mkdir -p /tmp/glob_logs
pids=()
for c in $(seq 1 $K); do
  Rscript c1_shapley_global_chunk.R "$c" "$K" > "/tmp/glob_logs/chunk_${c}.log" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} global chunks: ${pids[*]}"
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "ALL GLOBAL CHUNKS DONE (failures=$fail)"
ls out_files/Chapter3/nc_shapley2x_global | wc -l
ls out_files/Chapter3/tables/shapley_parts_global/ | wc -l
