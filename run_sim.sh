#!/bin/bash
set -euo pipefail
source ~/venvs/evo_sims/bin/activate


python run_evolution.py \
  --config-tsv simulation_config.tsv \
  --start-fasta start.fasta \
  --orf-file orfs.tsv \
  --reference-rscu reference.tsv \
  --penalty-regions penalty_regions.tsv \
  --output-dir output
