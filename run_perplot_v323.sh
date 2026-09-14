#!/usr/bin/env bash
# v3.2.3 per-plot sensitivity (ΔLAI/ΔHmax/ΔfCover + LAD real-vs-uniform) over the
# 400 cLHS plots, 12 chunks parallel. 8 sims/plot, resumable. ~90 min.
set -u
cd /home/corroyez/Documents/z_Example_rmusica_31012025
K=12
mkdir -p /tmp/perplot_v323_logs
pids=()
for c in $(seq 1 $K); do
  Rscript c1_sensitivity_perplot_chunk_v323.R "$c" "$K" > "/tmp/perplot_v323_logs/chunk_${c}.log" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} v323 chunks: ${pids[*]}"
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "ALL V323 PERPLOT CHUNKS DONE (failures=$fail)"
echo "nc: $(ls out_files/Chapter3/nc_sensitivity_perplot_v323 2>/dev/null | wc -l)"
echo "rows: $(cat out_files/Chapter3/tables/sensitivity_perplot_v323/part_*.csv 2>/dev/null | grep -vc '^pid')"
