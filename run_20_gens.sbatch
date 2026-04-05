#!/bin/bash
#SBATCH --job-name=evolve20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Safer bash behavior:
# -e        stop on command failure
# -u        stop on undefined variables
# pipefail  fail if any command in a pipeline fails
set -euo pipefail

# Required inputs
START_FASTA=${START_FASTA:?need START_FASTA}
ORF_FILE=${ORF_FILE:?need ORF_FILE}
REFERENCE_RSCU=${REFERENCE_RSCU:?need REFERENCE_RSCU}
CONFIG_TSV=${CONFIG_TSV:?need CONFIG_TSV}
PENALTY_REGIONS=${PENALTY_REGIONS:?need PENALTY_REGIONS}
OUTDIR=${OUTDIR:?need OUTDIR}

mkdir -p "$OUTDIR" logs

# Activate Python environment here if needed
source ~/venvs/evo_sims/bin/activate

python run_evolution.py \
  --config-tsv "$CONFIG_TSV" \
  --start-fasta "$START_FASTA" \
  --orf-file "$ORF_FILE" \
  --reference-rscu "$REFERENCE_RSCU" \
  --penalty-regions "PENALTY_REGIONS" \
  --output-dir "$OUTDIR"
