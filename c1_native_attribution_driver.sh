#!/bin/bash
# Native-20 m attribution: 400 plots x 8 sims = 3200 MuSICA runs, split into K chunks
# run in parallel. Resumable (run_musica_one skips existing nc). Detached.
set -u
K=10
LOGDIR=out_files/Chapter1; mkdir -p "$LOGDIR"
LOG="$LOGDIR/native20_attrib.log"; : > "$LOG"
echo "START $(date) K=$K" >> "$LOG"
seq 1 $K | xargs -P $K -I{} bash -c '
  Rscript c1_sensitivity_perplot_chunk_native20.R {} '"$K"' > '"$LOGDIR"'/native20_attrib_chunk_{}.log 2>&1
  echo "chunk {} done $(date)" >> "'"$LOG"'"
'
echo "ALL_DONE $(date)" >> "$LOG"
