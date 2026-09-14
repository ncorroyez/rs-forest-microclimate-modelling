#!/usr/bin/env bash
# iter-mode per-plot sensitivity (real MERRA-2 PBLH h_sbl), 400 cLHS, 12 chunks. 8 sims/plot, resumable.
set -u; cd /home/corroyez/Documents/z_Example_rmusica_31012025; K=12
mkdir -p /tmp/perplot_iter_logs; pids=()
for c in $(seq 1 $K); do Rscript c1_sensitivity_perplot_chunk_iter.R "$c" "$K" > "/tmp/perplot_iter_logs/chunk_${c}.log" 2>&1 & pids+=($!); done
echo "launched ${#pids[@]} iter chunks: ${pids[*]}"; fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "ALL ITER CHUNKS DONE (failures=$fail)"
echo "nc: $(ls out_files/Chapter3/nc_sensitivity_perplot_iter 2>/dev/null|wc -l) | rows: $(cat out_files/Chapter3/tables/sensitivity_perplot_iter/part_*.csv 2>/dev/null|grep -vc '^pid')"
