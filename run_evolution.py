import argparse
import csv
import math
import os
import random
from collections import Counter
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from evolve_worker import (
    get_mutation_probabilities,
    load_orf_regions,
    simulate_parent_batch,
)

# User-specified amino acid order
AA_ORDER = (
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "Stop"
)

# User-specified codon grouping and ordering
CODON_GROUPS_ORDERED = [
    ["GCT", "GCC", "GCA", "GCG"],                         # A
    ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],          # R
    ["AAT", "AAC"],                                       # N
    ["GAT", "GAC"],                                       # D
    ["TGT", "TGC"],                                       # C
    ["CAA", "CAG"],                                       # Q
    ["GAA", "GAG"],                                       # E
    ["GGT", "GGC", "GGA", "GGG"],                         # G
    ["CAT", "CAC"],                                       # H
    ["ATT", "ATC", "ATA"],                                # I
    ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],          # L
    ["AAA", "AAG"],                                       # K
    ["ATG"],                                              # M
    ["TTT", "TTC"],                                       # F
    ["CCT", "CCC", "CCA", "CCG"],                         # P
    ["TCA", "TCC", "TCT", "TCG", "AGT", "AGC"],          # S
    ["ACA", "ACT", "ACG", "ACC"],                         # T
    ["TGG"],                                              # W
    ["TAC", "TAT"],                                       # Y
    ["GTA", "GTC", "GTT", "GTG"],                         # V
    ["TGA", "TAG", "TAA"]                                 # Stop
]

AA_TO_CODONS = {
    aa: codons
    for aa, codons in zip(AA_ORDER, CODON_GROUPS_ORDERED)
}

CODON_TO_AA = {}
for aa, codons in AA_TO_CODONS.items():
    for codon in codons:
        CODON_TO_AA[codon] = aa

# Flat codon order, including stop codons
CODON_ORDER = [
    codon
    for group in CODON_GROUPS_ORDERED
    for codon in group
]

# Sense codons only for RSCU / reference comparison
SENSE_AA_TO_CODONS = {
    aa: codons
    for aa, codons in AA_TO_CODONS.items()
    if aa != "Stop"
}

SENSE_CODON_ORDER = [
    codon
    for aa in AA_ORDER
    if aa != "Stop"
    for codon in AA_TO_CODONS[aa]
]


def load_config_tsv(path):
    """
    Load a simple key-value TSV with columns:
        parameter   value

    Values are automatically cast to int or float when possible.
    """
    config = {}

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"parameter", "value"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(f"Config TSV must contain columns: {sorted(required)}")

        for row in reader:
            key = row["parameter"]
            value = row["value"]

            if value.isdigit():
                value = int(value)
            else:
                try:
                    value = float(value)
                except ValueError:
                    pass

            config[key] = value

    return config


def load_initial_parent(start_fasta: str):
    """
    Load the initial sequence from FASTA and wrap it as the first parent.
    """
    record = next(SeqIO.parse(start_fasta, "fasta"))
    return [{
        "parent_id": record.id,
        "sequence": str(record.seq),
        "stop_codon_count": 0,
        "is_viable": True,
        "stop_codon": "",
        "codon_start_1based": "",
        "codon_index_global_1based": "",
        "first_stop_orf_name": "",
        "codon_index_within_orf_1based": "",
        "mutated_nt_index_1based": "",
        "reference_used": "",
        "similarity_index": "",
        "reproductive_score": ""
    }]


def parent_can_reproduce(parent, tolerated_stop_codons: int) -> bool:
    """
    Decide whether a parent can produce descendants.

    A parent reproduces if its stop-codon count is <= threshold.
    """
    return parent["stop_codon_count"] <= tolerated_stop_codons


def extract_concatenated_orf_sequence(sequence: str, orf_regions):
    """
    Extract ORFs from the sequence and concatenate them in ORF-file order.

    This concatenated coding sequence is used for codon counting and RSCU.
    """
    parts = []
    for orf in orf_regions:
        start0 = orf["start_1based"] - 1
        end0_exclusive = orf["end_1based"]
        parts.append(sequence[start0:end0_exclusive])
    return "".join(parts)


def count_codons_in_sequence(coding_sequence: str):
    """
    Count codons in a coding sequence whose length must be divisible by 3.
    """
    if len(coding_sequence) % 3 != 0:
        raise ValueError("Concatenated ORF sequence length is not divisible by 3.")

    codon_counts = Counter()
    for i in range(0, len(coding_sequence), 3):
        codon = coding_sequence[i:i+3]
        codon_counts[codon] += 1
    return codon_counts


def compute_rscu_vector(sequence: str, orf_regions):
    """
    Compute an RSCU vector in the fixed codon order requested by the user.

    Stop codons are excluded from the RSCU vector.

    For each synonymous codon group:
        RSCU(codon) = observed_count / expected_count_if_uniform
    """
    coding_sequence = extract_concatenated_orf_sequence(sequence, orf_regions)
    codon_counts = count_codons_in_sequence(coding_sequence)

    rscu = {}

    for aa, codons in SENSE_AA_TO_CODONS.items():
        total = sum(codon_counts.get(codon, 0) for codon in codons)
        n_syn = len(codons)

        if total == 0:
            # If no codons from this family appear, assign 0 to all
            for codon in codons:
                rscu[codon] = 0.0
        else:
            expected = total / n_syn
            for codon in codons:
                rscu[codon] = codon_counts.get(codon, 0) / expected

    return [rscu[codon] for codon in SENSE_CODON_ORDER]


def cosine_distance(vec1, vec2):
    """
    Compute cosine distance between two vectors.

    In this project, cosine distance is the Similarity Index:
      lower = more similar
      higher = less similar
    """
    if len(vec1) != len(vec2):
        raise ValueError("Vectors must be same length.")

    dot = sum(a * b for a, b in zip(vec1, vec2))
    norm1 = math.sqrt(sum(a * a for a in vec1))
    norm2 = math.sqrt(sum(b * b for b in vec2))

    if norm1 == 0.0 or norm2 == 0.0:
        return 1.0

    cosine_similarity = dot / (norm1 * norm2)
    return 1.0 - cosine_similarity


def load_reference_rscu_tsv(reference_tsv: str):
    """
    Load reference RSCU vectors from a TSV.

    Required columns:
      reference_name
      one column for every codon in CODON_ORDER

    Stop-codon columns may be present in the file, but only sense codons
    are used in the comparison vectors.
    """
    references = {}

    with open(reference_tsv, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"reference_name"} | set(CODON_ORDER)
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Reference TSV missing required columns: {sorted(missing)}")

        for row in reader:
            name = row["reference_name"]
            references[name] = {
                "sense_vector": [float(row[codon]) for codon in SENSE_CODON_ORDER]
            }

    if len(references) < 2:
        raise ValueError("Reference TSV must contain at least two reference rows.")

    return references


def get_reference_name_for_generation(generation: int, cycle_n: int, ref1_name: str, ref2_name: str):
    """
    Determine which reference is active for a generation.

    Rule:
    - first n generations in each cycle use reference 1
    - generation n+1 uses reference 2
    - then the cycle repeats

    Example if n = 4:
      1,2,3,4 -> ref1
      5       -> ref2
      6,7,8,9 -> ref1
      10      -> ref2
    """
    cycle_len = cycle_n + 1
    pos = (generation - 1) % cycle_len
    if pos < cycle_n:
        return ref1_name
    return ref2_name


def compute_scalings(ranks, adaptation_scaling, mode):
    """
    Convert parent ranks into scaling factors.

    rank 0 = best parent
    rank N-1 = worst parent

    Modes:
      linear
      exponential
      exponential_normalized
    """
    N = len(ranks)

    if N <= 1:
        return [1.0]

    mean_rank = (N - 1) / 2.0
    normalized_ranks = [
        (rank - mean_rank) / mean_rank
        for rank in ranks
    ]

    if mode == "linear":
        scalings = [
            max(0.0, 1.0 - adaptation_scaling * r)
            for r in normalized_ranks
        ]

    elif mode == "exponential":
        scalings = [
            math.exp(-adaptation_scaling * r)
            for r in normalized_ranks
        ]

    elif mode == "exponential_normalized":
        raw = [
            math.exp(-adaptation_scaling * r)
            for r in normalized_ranks
        ]
        mean_val = sum(raw) / len(raw)
        scalings = [x / mean_val for x in raw]

    else:
        raise ValueError(f"Unknown adaptation mode: {mode}")

    return scalings


def write_parent_set(parents, outdir):
    """
    Write the sampled parent set for a generation.

    Outputs:
    - one FASTA per parent
    - parent manifest
    - parent metadata TSV
    """
    os.makedirs(outdir, exist_ok=True)

    manifest_path = os.path.join(outdir, "parents_manifest.txt")
    metadata_path = os.path.join(outdir, "parent_metadata.tsv")

    with open(manifest_path, "w") as manifest, open(metadata_path, "w", newline="") as meta:
        writer = csv.writer(meta, delimiter="\t")
        writer.writerow([
            "parent_id",
            "stop_codon_count",
            "is_viable",
            "stop_codon",
            "codon_start_1based",
            "codon_index_global_1based",
            "first_stop_orf_name",
            "codon_index_within_orf_1based",
            "mutated_nt_index_1based",
            "reference_used",
            "similarity_index",
            "reproductive_score",
            "sequence"
        ])

        for i, parent in enumerate(parents):
            fasta_id = f"parent_{i:03d}"
            fasta_path = os.path.join(outdir, f"{fasta_id}.fasta")

            rec = SeqRecord(Seq(parent["sequence"]), id=fasta_id, description="")
            with open(fasta_path, "w") as handle:
                SeqIO.write(rec, handle, "fasta")

            manifest.write(fasta_path + "\n")

            writer.writerow([
                parent["parent_id"],
                parent["stop_codon_count"],
                int(parent["is_viable"]),
                parent["stop_codon"],
                parent["codon_start_1based"],
                parent["codon_index_global_1based"],
                parent["first_stop_orf_name"],
                parent["codon_index_within_orf_1based"],
                parent["mutated_nt_index_1based"],
                parent["reference_used"],
                parent["similarity_index"],
                parent["reproductive_score"],
                parent["sequence"]
            ])


def write_checkpoint_record(record_dir, generation, variants):
    """
    Write full checkpoint outputs for selected generations.
    """
    os.makedirs(record_dir, exist_ok=True)

    viable_fasta = os.path.join(record_dir, f"gen_{generation:03d}_viable_variants.fasta")
    stop_fasta = os.path.join(record_dir, f"gen_{generation:03d}_stop_variants.fasta")
    variant_log = os.path.join(record_dir, f"gen_{generation:03d}_all_variants.tsv")

    with open(viable_fasta, "w") as viable_out, open(stop_fasta, "w") as stop_out, open(variant_log, "w", newline="") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t")
        writer.writerow([
            "simulation_index",
            "source_parent_id",
            "stop_codon_count",
            "is_viable",
            "stop_codon",
            "codon_start_1based",
            "codon_index_global_1based",
            "first_stop_orf_name",
            "codon_index_within_orf_1based",
            "mutated_nt_index_1based",
            "sequence"
        ])

        viable_i = 1
        stop_i = 1

        for variant in variants:
            writer.writerow([
                variant["simulation_index"],
                variant["source_parent_id"],
                variant["stop_codon_count"],
                int(variant["is_viable"]),
                variant["stop_codon"],
                variant["codon_start_1based"],
                variant["codon_index_global_1based"],
                variant["first_stop_orf_name"],
                variant["codon_index_within_orf_1based"],
                variant["mutated_nt_index_1based"],
                variant["sequence"]
            ])

            if variant["stop_codon_count"] == 0:
                viable_out.write(f">gen{generation:03d}_viable_{viable_i}\n{variant['sequence']}\n")
                viable_i += 1
            else:
                stop_out.write(f">gen{generation:03d}_stop_{stop_i}\n{variant['sequence']}\n")
                stop_i += 1


def append_generation_summary(summary_path, generation, reference_used, sampled_parent_count,
                              reproductive_parent_count, total_requested_progeny,
                              total_realized_variants, fully_viable_variant_count,
                              stop_variant_count, taa_count, tag_count, tga_count,
                              mean_similarity_index_ref1, mean_similarity_index_ref2,
                              checkpoint_written):
    """
    Append one row of summary metrics for the current generation.
    """
    write_header = not os.path.exists(summary_path)

    with open(summary_path, "a", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        if write_header:
            writer.writerow([
                "generation",
                "reference_used",
                "sampled_parent_count",
                "reproductive_parent_count",
                "total_requested_progeny",
                "total_realized_variants",
                "fully_viable_variant_count",
                "stop_variant_count",
                "TAA_count",
                "TAG_count",
                "TGA_count",
                "mean_similarity_index_ref1",
                "mean_similarity_index_ref2",
                "checkpoint_written"
            ])

        writer.writerow([
            generation,
            reference_used,
            sampled_parent_count,
            reproductive_parent_count,
            total_requested_progeny,
            total_realized_variants,
            fully_viable_variant_count,
            stop_variant_count,
            taa_count,
            tag_count,
            tga_count,
            mean_similarity_index_ref1,
            mean_similarity_index_ref2,
            int(checkpoint_written)
        ])


def main():
    parser = argparse.ArgumentParser()

    # File inputs
    parser.add_argument("--config-tsv", required=True)
    parser.add_argument("--start-fasta", required=True)
    parser.add_argument("--orf-file", required=True)
    parser.add_argument("--reference-rscu", required=True)
    parser.add_argument("--output-dir", required=True)

    args = parser.parse_args()

    # Load config and required resources
    cfg = load_config_tsv(args.config_tsv)
    orf_regions = load_orf_regions(args.orf_file)
    mutation_probabilities = get_mutation_probabilities(int(cfg["mutation_set"]))
    reference_vectors = load_reference_rscu_tsv(args.reference_rscu)

    # Pull parameters from config
    mutation_rate = float(cfg["mutation_rate"])
    num_generations = int(cfg["num_generations"])
    record_interval = int(cfg["record_interval"])
    sample_size = int(cfg["sample_size"])
    progeny_scale = float(cfg["progeny_scale"])
    tolerated_stop_codons = int(cfg["tolerated_stop_codons"])
    cpus = int(cfg["cpus"])
    seed = int(cfg["seed"])

    reference1_name = cfg["reference1_name"]
    reference2_name = cfg["reference2_name"]
    reference_cycle_n = int(cfg["reference_cycle_n"])

    adaptation_scaling_ref1 = float(cfg["adaptation_scaling_ref1"])
    adaptation_scaling_ref2 = float(cfg["adaptation_scaling_ref2"])
    adaptation_mode = cfg["adaptation_mode"]

    max_progeny_per_parent = int(cfg["max_progeny_per_parent"])

    # Validate reference names
    if reference1_name not in reference_vectors:
        raise ValueError(f"Reference not found in reference TSV: {reference1_name}")
    if reference2_name not in reference_vectors:
        raise ValueError(f"Reference not found in reference TSV: {reference2_name}")

    # Random generator for parent sampling between generations
    rng = random.Random(seed)

    os.makedirs(args.output_dir, exist_ok=True)

    # Load initial parent and write generation-1 parent files
    current_parents = load_initial_parent(args.start_fasta)
    write_parent_set(current_parents, os.path.join(args.output_dir, "gen_001_parents"))

    summary_path = os.path.join(args.output_dir, "generation_summaries.tsv")

    # Main evolutionary loop
    for generation in range(1, num_generations + 1):
        # Determine which reference is active this generation
        reference_name = get_reference_name_for_generation(
            generation=generation,
            cycle_n=reference_cycle_n,
            ref1_name=reference1_name,
            ref2_name=reference2_name
        )

        ref1_vector = reference_vectors[reference1_name]["sense_vector"]
        ref2_vector = reference_vectors[reference2_name]["sense_vector"]
        active_reference_vector = reference_vectors[reference_name]["sense_vector"]

        # Reproductive parents are those whose stop count is <= threshold
        reproductive_parents = [
            p for p in current_parents
            if parent_can_reproduce(p, tolerated_stop_codons)
        ]

        print(f"Generation {generation}: sampled parents = {len(current_parents)}")
        print(f"Generation {generation}: reproductive parents = {len(reproductive_parents)}")
        print(f"Generation {generation}: active reference = {reference_name}")

        if not reproductive_parents:
            raise RuntimeError(
                f"No parents available to produce generation {generation} "
                f"under current stop-codon threshold."
            )

        # Compute Similarity Index to both references for logging,
        # and to the active reference for ranking.
        parent_scores = []
        ref1_distances = []
        ref2_distances = []

        for parent in reproductive_parents:
            parent_rscu = compute_rscu_vector(parent["sequence"], orf_regions)

            dist_ref1 = cosine_distance(parent_rscu, ref1_vector)
            dist_ref2 = cosine_distance(parent_rscu, ref2_vector)

            ref1_distances.append(dist_ref1)
            ref2_distances.append(dist_ref2)

            # Active reference determines selection
            active_dist = cosine_distance(parent_rscu, active_reference_vector)

            parent["reference_used"] = reference_name
            parent["similarity_index"] = active_dist

            parent_scores.append((parent, active_dist))

        mean_similarity_index_ref1 = sum(ref1_distances) / len(ref1_distances)
        mean_similarity_index_ref2 = sum(ref2_distances) / len(ref2_distances)

        # Lower cosine distance means more similar, so sort ascending
        parent_scores.sort(key=lambda x: x[1])

        # Select adaptation scaling based on active reference
        if reference_name == reference1_name:
            adaptation_scaling = adaptation_scaling_ref1
        else:
            adaptation_scaling = adaptation_scaling_ref2

        ranks = list(range(len(parent_scores)))
        scalings = compute_scalings(
            ranks=ranks,
            adaptation_scaling=adaptation_scaling,
            mode=adaptation_mode
        )

        # Convert scaling factors into progeny counts
        progeny_counts = []
        for (parent, dist), scaling in zip(parent_scores, scalings):
            parent["reproductive_score"] = scaling

            progeny = int(round(progeny_scale * scaling))
            progeny = min(progeny, max_progeny_per_parent)

            progeny_counts.append((parent, progeny))

        total_requested_progeny = sum(n for _, n in progeny_counts)

        if total_requested_progeny == 0:
            raise RuntimeError(
                f"All reproductive parents received zero progeny in generation {generation}."
            )

        # Build parallel simulation tasks
        task_args = []
        for parent_idx, (parent, n_sims) in enumerate(progeny_counts):
            if n_sims <= 0:
                continue

            task_args.append((
                parent["parent_id"],
                parent["sequence"],
                n_sims,
                mutation_rate,
                mutation_probabilities,
                orf_regions,
                seed + generation * 100000 + parent_idx
            ))

        # Simulate all descendants in parallel
        all_variants = []
        with Pool(processes=cpus) as pool:
            for variant_list in pool.map(simulate_parent_batch, task_args):
                all_variants.extend(variant_list)

        if not all_variants:
            raise RuntimeError(f"No variants produced in generation {generation}.")

        # Summaries of generated variants
        fully_viable_variants = [v for v in all_variants if v["stop_codon_count"] == 0]
        stop_variants = [v for v in all_variants if v["stop_codon_count"] > 0]

        stop_counter = Counter(v["stop_codon"] for v in stop_variants if v["stop_codon"])
        taa_count = stop_counter.get("TAA", 0)
        tag_count = stop_counter.get("TAG", 0)
        tga_count = stop_counter.get("TGA", 0)

        # Write checkpoint outputs every record_interval generations
        checkpoint_written = (generation % record_interval == 0)
        if checkpoint_written:
            record_dir = os.path.join(args.output_dir, f"gen_{generation:03d}_record")
            write_checkpoint_record(record_dir, generation, all_variants)

        # Sample next generation parents from ALL variants, including stop variants
        if generation < num_generations:
            sampled_n = min(sample_size, len(all_variants))
            sampled_variants = rng.sample(all_variants, sampled_n)

            next_parents = []
            for i, variant in enumerate(sampled_variants):
                next_parents.append({
                    "parent_id": f"gen{generation + 1:03d}_parent_{i:03d}",
                    "sequence": variant["sequence"],
                    "stop_codon_count": variant["stop_codon_count"],
                    "is_viable": variant["is_viable"],
                    "stop_codon": variant["stop_codon"],
                    "codon_start_1based": variant["codon_start_1based"],
                    "codon_index_global_1based": variant["codon_index_global_1based"],
                    "first_stop_orf_name": variant["first_stop_orf_name"],
                    "codon_index_within_orf_1based": variant["codon_index_within_orf_1based"],
                    "mutated_nt_index_1based": variant["mutated_nt_index_1based"],
                    "reference_used": "",
                    "similarity_index": "",
                    "reproductive_score": ""
                })

            next_parent_dir = os.path.join(args.output_dir, f"gen_{generation + 1:03d}_parents")
            write_parent_set(next_parents, next_parent_dir)
        else:
            next_parents = []

        # Write one generation summary row
        append_generation_summary(
            summary_path=summary_path,
            generation=generation,
            reference_used=reference_name,
            sampled_parent_count=len(current_parents),
            reproductive_parent_count=len(reproductive_parents),
            total_requested_progeny=total_requested_progeny,
            total_realized_variants=len(all_variants),
            fully_viable_variant_count=len(fully_viable_variants),
            stop_variant_count=len(stop_variants),
            taa_count=taa_count,
            tag_count=tag_count,
            tga_count=tga_count,
            mean_similarity_index_ref1=mean_similarity_index_ref1,
            mean_similarity_index_ref2=mean_similarity_index_ref2,
            checkpoint_written=checkpoint_written
        )

        # Move forward one generation
        current_parents = next_parents

    print(f"Done. Summary written to {summary_path}")


if __name__ == "__main__":
    main()
