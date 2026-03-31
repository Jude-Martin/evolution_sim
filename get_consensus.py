#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def compute_consensus_and_stats(sequences):
    """
    Build a consensus sequence and per-position statistics.

    Returns:
        consensus_sequence (str)
        stats_rows (list of dicts)
    To run use     
    python consensus_with_stats.py \
  --input-fasta input.fasta \
  --output-consensus-fasta consensus.fasta \
  --output-stats-tsv consensus_stats.tsv
    
    """
    if not sequences:
        raise ValueError("No sequences found.")

    seq_length = len(sequences[0])

    for seq in sequences:
        if len(seq) != seq_length:
            raise ValueError("All sequences must have the same length.")

    consensus = []
    stats_rows = []

    for i in range(seq_length):
        column = [seq[i].upper() for seq in sequences]
        counts = Counter(column)

        total = sum(counts.values())

        # Ensure all bases are represented
        base_counts = {b: counts.get(b, 0) for b in ["A", "C", "G", "T", "N"]}

        # Determine consensus (handle ties)
        most_common = counts.most_common()
        if len(most_common) == 1 or most_common[0][1] > most_common[1][1]:
            consensus_base = most_common[0][0]
        else:
            consensus_base = "N"

        consensus.append(consensus_base)

        stats_rows.append({
            "position_1based": i + 1,
            "consensus_base": consensus_base,
            "A_count": base_counts["A"],
            "C_count": base_counts["C"],
            "G_count": base_counts["G"],
            "T_count": base_counts["T"],
            "N_count": base_counts["N"],
            "A_freq": base_counts["A"] / total,
            "C_freq": base_counts["C"] / total,
            "G_freq": base_counts["G"] / total,
            "T_freq": base_counts["T"] / total,
            "N_freq": base_counts["N"] / total,
            "coverage": total
        })

    return "".join(consensus), stats_rows


def write_stats_tsv(stats_rows, output_path):
    """
    Write per-position statistics to TSV.
    """
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "position_1based",
                "consensus_base",
                "A_count", "C_count", "G_count", "T_count", "N_count",
                "A_freq", "C_freq", "G_freq", "T_freq", "N_freq",
                "coverage"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(stats_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Generate consensus sequence and per-position nucleotide statistics from FASTA."
    )
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-consensus-fasta", required=True)
    parser.add_argument("--output-stats-tsv", required=True)
    parser.add_argument("--consensus-id", default="consensus")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    sequences = [str(record.seq) for record in records]

    consensus_seq, stats_rows = compute_consensus_and_stats(sequences)

    consensus_record = SeqRecord(
        Seq(consensus_seq),
        id=args.consensus_id,
        description=f"consensus_from={args.input_fasta} n_sequences={len(sequences)}"
    )

    with open(args.output_consensus_fasta, "w") as handle:
        SeqIO.write(consensus_record, handle, "fasta")

    write_stats_tsv(stats_rows, args.output_stats_tsv)

    print(f"Wrote consensus FASTA: {args.output_consensus_fasta}")
    print(f"Wrote per-position stats TSV: {args.output_stats_tsv}")


if __name__ == "__main__":
    main()
