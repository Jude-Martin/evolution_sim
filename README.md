Simple simulation of a 2 host system with bottlenecks in VEEV fidelity/mutator variants and stop mutants.


Bottlenecks will influence drift.  Selection is determined by RSCU values (comparing to reference file).
Parameters can be set in the orfs.tsv and simulation_config.tsv files.

Sequence starts from a single entry fasta file. 
Various parameters can be changed:

Number of generations to run
Number of replication events to simulate per parent
Number of parents to sample per generation (simulating drft/bottlenecks)
Influence exertd by each host
Number of generations in host 1 before switch to host 2
position of ORFs in starting sequence.
Number of stop codons tolerated, should be at least equal to number of ORFs, can be higher

Mutational frequencies are specified in evolve_worker.py 

To run simulation:
sbatch --export=ALL,START_FASTA=start.fasta,ORF_FILE=orfs.tsv,REFERENCE_RSCU=reference.tsv,CONFIG_TSV=simulation_config.tsv,OUTDIR=output run_20gens.sbatch

To plot some metrics:

python plot_evolution_metrics.py \
  --summary-tsv output/generation_summaries.tsv \
  --output-dir output/plots

Currently 3 modes can be used for positive selection, exponential_normalized, linear, exponential
based on scoring of similarity of each parent to reference and scaling by a specified factor

mean_rank is calculated, normalised converted to a weighting factor and used to calculate # of progeny
	​


To check orfs prior to running:

python check_orfs.py   --start-fasta start.fasta --orf-file orfs.tsv  --output-fasta extracted_orfs.fasta

To plot hotspots for stop codon sites

python plot_stop_hotspots.py \
  --input-tsv output/gen_005_record/gen_005_novel_stop_events.tsv \
  --label gen005 \
  --input-tsv output/gen_010_record/gen_010_novel_stop_events.tsv \
  --label gen010 \
  --sequence-length 2400 \
  --window-width 50 \
  --shade-regions shade_regions.tsv \
  --output-dir hotspot_plots
