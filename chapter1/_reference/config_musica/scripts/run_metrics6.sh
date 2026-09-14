#!/usr/bin/env bash
set -u; cd /home/corroyez/Documents/z_Example_rmusica_31012025; K=12
mkdir -p /tmp/metrics6_logs
for VER in v320 iter; do
  pids=()
  for c in $(seq 1 $K); do Rscript c1_metrics6_perplot.R "$c" "$K" "$VER" > "/tmp/metrics6_logs/${VER}_${c}.log" 2>&1 & pids+=($!); done
  for p in "${pids[@]}"; do wait "$p"; done
  echo "$VER DONE: rows=$(cat out_files/Chapter3/tables/metrics6_${VER}/part_*.csv 2>/dev/null|grep -vc '^pid')"
done
echo "ALL METRICS6 DONE"
