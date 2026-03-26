Simple simulation of a 2 host system with bottlenecks in VEEV fidelity/mutator variants and stop mutants.


Bottlenecks will influence drift.  Selection is determined by RSCU values (comparing to reference file).
Parameters can be set in the orfs.tsv and simulation_config.tsv files. 



To run simulation:
sbatch --export=ALL,START_FASTA=start.fasta,ORF_FILE=orfs.tsv,REFERENCE_RSCU=reference.tsv,CONFIG_TSV=simulation_config.tsv,OUTDIR=output run_20gens.sbatch

To plot some metrics:

python plot_evolution_metrics.py \
  --summary-tsv output/generation_summaries.tsv \
  --output-dir output/plots
