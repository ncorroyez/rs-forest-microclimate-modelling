#!/usr/bin/env bash
# Relance les 12 chunks per-cluster Shapley (legacy binary). Reprenable: les nc
# existants sont relus par run_musica_one (file.exists -> return). Tourne les 12
# en parallele, attend tous, log par chunk.
set -u
cd /home/corroyez/Documents/z_Example_rmusica_31012025
K=12
mkdir -p /tmp/pcl_logs
pids=()
for c in $(seq 1 $K); do
  Rscript c3_shapley_chunk.R "$c" "$K" > "/tmp/pcl_logs/chunk_${c}.log" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} chunks: ${pids[*]}"
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "ALL CHUNKS DONE (failures=$fail)"
ls out_files/Chapter3/nc_shapley2x | wc -l
ls out_files/Chapter3/tables/shapley_parts/ | wc -l
