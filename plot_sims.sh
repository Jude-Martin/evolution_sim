#!/bin/bash
set -euo pipefail
source ~/venvs/evo_sims/bin/activate
OUTDIR=combined_plots

python plot_evolution_metrics.py \
  --summary-tsv run1/output/generation_summaries.tsv \
  --label run1 \
  --summary-tsv run2/output/generation_summaries.tsv \
  --label run2 \
  --summary-tsv run3/output/generation_summaries.tsv \
  --label run3 \
  --output-dir "$OUTDIR" \
  --plots similarity dnds viable_parent_count novel_stop_count \
  --ymax-similarity 0.5 \
  --ymax-dnds 3.0 \
  --ymax-viable-parent-count 100 \
  --ymax-novel-stop-count 100000
