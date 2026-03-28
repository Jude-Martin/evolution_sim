#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def load_orf_regions(orffile: str):
    """
    Load ORFs from a TSV with columns:
        orf_name    start_1based    end_1based

    Coordinates are 1-based inclusive.
    Each ORF length must be divisible by 3.
    """
    orfs = []

    with open(orffile, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"orf_name", "start_1based", "end_1based"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(
                f"ORF TSV must contain columns: {sorted(required)}"
            )

        for row in reader:
            orf_name = row["orf_name"]
            start = int(row["start_1based"])
            end = int(row["end_1based"])

            if start < 1 or end < start:
                raise ValueError(
                    f"Invalid coordinates for {orf_name}: {start}-{end}"
                )

            length = end - start + 1
            if length % 3 != 0:
                raise ValueError(
                    f"ORF {orf_name} length is not divisible by 3: "
                    f"{length} nt ({start}-{end})"
                )

            orfs.append({
                "orf_name": orf_name,
                "start_1based": start,
                "end_1based": end,
            })

    if not orfs:
        raise ValueError("No ORFs found in ORF TSV.")

    return orfs


def extract_orfs_from_sequence(sequence: str, sequence_id: str, orfs):
    """
    Extract each ORF as a separate SeqRecord.
    """
    records = []
    seq_len = len(sequence)

    for orf in orfs:
        start_1based = orf["start_1based"]
        end_1based = orf["end_1based"]

        if end_1based > seq_len:
            raise ValueError(
                f"ORF {orf['orf_name']} ends beyond sequence length: "
                f"{end_1based} > {seq_len}"
            )

        start0 = start_1based - 1
        end0_exclusive = end_1based
        orf_seq = sequence[start0:end0_exclusive]

        record_id = f"{sequence_id}|{orf['orf_name']}"
        description = (
            f"orf_name={orf['orf_name']} "
            f"start_1based={start_1based} "
            f"end_1based={end_1based} "
            f"length={len(orf_seq)}"
        )

        records.append(
            SeqRecord(
                Seq(orf_seq),
                id=record_id,
                description=description,
            )
        )

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Extract ORFs from a starting FASTA and write each ORF as a separate FASTA record."
    )
    parser.add_argument(
        "--start-fasta",
        required=True,
        help="Input FASTA containing the starting sequence",
    )
    parser.add_argument(
        "--orf-file",
        required=True,
        help="TSV file with columns: orf_name, start_1based, end_1based",
    )
    parser.add_argument(
        "--output-fasta",
        required=True,
        help="Output FASTA file containing one record per ORF",
    )
    args = parser.parse_args()

    # Read input sequence
    record = next(SeqIO.parse(args.start_fasta, "fasta"))
    sequence = str(record.seq)
    sequence_id = record.id

    # Load ORFs
    orfs = load_orf_regions(args.orf_file)

    # Extract ORFs
    orf_records = extract_orfs_from_sequence(sequence, sequence_id, orfs)

    # Write output
    output_path = Path(args.output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as handle:
        SeqIO.write(orf_records, handle, "fasta")

    print(f"Wrote {len(orf_records)} ORFs to {output_path}")


if __name__ == "__main__":
    main()
