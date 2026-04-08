import csv
import random

# Standard DNA stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Mutation probability matrices.
mutation_probabilities_set1_68U201 = {
    'A': [0.0406342 / 0.22313, 0.1569820 / 0.22313, 0.0255137 / 0.22313],
    'C': [0.0715498 / 0.258425, 0.1093084 / 0.258425, 0.0775670 / 0.258425],
    'G': [0.069933 / 0.138767, 0.0232894 / 0.138767, 0.0455444 / 0.138767],
    'T': [0.08778562 / 0.37967764, 0.08651281 / 0.37967764, 0.20537921 / 0.37967764]
}

mutation_probabilities_set2_G7R = {
    'A': [0.07799078 / 0.24349108, 0.13299544 / 0.24349108, 0.03250486 / 0.24349108],
    'C': [0.08857339 / 0.239979, 0.08100787 / 0.239979, 0.07039774 / 0.239979],
    'G': [0.07645003 / 0.19021658, 0.08419321 / 0.19021658, 0.02957334 / 0.19021658],
    'T': [0.07249201 / 0.32631333, 0.10876877 / 0.32631333, 0.14505255 / 0.32631333]
}

mutation_probabilities_set3_3X = {
    'A': [0.03321736 / 0.19424901, 0.1352778 / 0.19424901, 0.02575385 / 0.19424901],
    'C': [0.07334453 / 0.22282149, 0.08144008 / 0.22282149, 0.06803688 / 0.22282149],
    'G': [0.07667533 / 0.1873479, 0.09463397 / 0.1873479, 0.0160386 / 0.1873479],
    'T': [0.14717589 / 0.39558158, 0.10889809 / 0.39558158, 0.13950761 / 0.39558158]
}

mutation_probabilities_set4_4X = {
    'A': [0.02102178 / 0.17761079, 0.13136807 / 0.17761079, 0.02522094 / 0.17761079],
    'C': [0.06261927 / 0.20481857, 0.08215234 / 0.20481857, 0.06004696 / 0.20481857],
    'G': [0.13718076 / 0.24754165, 0.09510129 / 0.24754165, 0.0152596 / 0.24754165],
    'T': [0.13576939 / 0.370029, 0.10137649 / 0.370029, 0.13288312 / 0.370029]
}

mutation_probabilities_set5_neutral = {
    'A': [1 / 3, 1 / 3, 1 / 3],
    'C': [1 / 3, 1 / 3, 1 / 3],
    'G': [1 / 3, 1 / 3, 1 / 3],
    'T': [1 / 3, 1 / 3, 1 / 3]
}


def get_mutation_probabilities(set_id: int):
    if set_id == 1:
        return mutation_probabilities_set1_68U201
    if set_id == 2:
        return mutation_probabilities_set2_G7R
    if set_id == 3:
        return mutation_probabilities_set3_3X
    if set_id == 4:
        return mutation_probabilities_set4_4X
     if set_id == 5:
        return mutation_probabilities_set5_neutral
    
    raise ValueError(f"Invalid mutation set: {set_id}")


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
            raise ValueError(f"ORF TSV must contain columns: {sorted(required)}")

        for row in reader:
            start = int(row["start_1based"])
            end = int(row["end_1based"])

            if start < 1 or end < start:
                raise ValueError(f"Invalid ORF coordinates: {row}")

            length = end - start + 1
            if length % 3 != 0:
                raise ValueError(
                    f"ORF length is not divisible by 3: {row['orf_name']} "
                    f"{start}-{end} ({length} nt)"
                )

            orfs.append({
                "orf_name": row["orf_name"],
                "start_1based": start,
                "end_1based": end
            })

    if not orfs:
        raise ValueError("No ORFs found in ORF TSV.")

    return orfs


def mutate_nucleotide(nucleotide: str, mutation_probabilities: dict) -> str:
    mutation_outcomes = {
        'A': ['C', 'G', 'T'],
        'C': ['A', 'G', 'T'],
        'G': ['A', 'C', 'T'],
        'T': ['A', 'C', 'G']
    }

    return random.choices(
        mutation_outcomes[nucleotide],
        mutation_probabilities[nucleotide]
    )[0]


def analyze_orf_stops(sequence: str, orfs):
    """
    Analyze all stop codons inside ORFs.

    Returns a stop list with:
      stop_codon
      codon_start_1based
      codon_index_global_1based
      first_stop_orf_name
      codon_index_within_orf_1based
      is_terminal_orf_stop
    """
    stops = []

    for orf in orfs:
        start0 = orf["start_1based"] - 1
        end1 = orf["end_1based"]
        terminal_start0 = end1 - 3

        codon_idx_within_orf = 0
        for i in range(start0, end1 - 2, 3):
            codon_idx_within_orf += 1
            codon = sequence[i:i+3]

            if codon in STOP_CODONS:
                stops.append({
                    "stop_codon": codon,
                    "codon_start_1based": i + 1,
                    "codon_index_global_1based": (i // 3) + 1,
                    "first_stop_orf_name": orf["orf_name"],
                    "codon_index_within_orf_1based": codon_idx_within_orf,
                    "is_terminal_orf_stop": (i == terminal_start0),
                })

    return {
        "total_stop_codon_count": len(stops),
        "stop_list": stops,
    }


def build_stop_signature_set(stop_list):
    """
    Convert a stop list into a set of stop-position signatures:
      (orf_name, codon_index_within_orf_1based)

    Stop identity is ignored.
    """
    signatures = set()
    for stop in stop_list:
        signatures.add((
            stop["first_stop_orf_name"],
            stop["codon_index_within_orf_1based"],
        ))
    return signatures


def summarize_stops_against_expected(sequence: str, orfs, expected_stop_signatures):
    """
    Compare a sequence's stop positions against the expected stop-position set
    derived from the starting FASTA.

    Expected signatures are position-only:
      (orf_name, codon_index_within_orf_1based)

    Returns:
      total_stop_codon_count
      novel_stop_codon_count
      missing_expected_stop_count
      stop_list
    """
    stop_info = analyze_orf_stops(sequence, orfs)
    stop_list = stop_info["stop_list"]

    observed_signatures = build_stop_signature_set(stop_list)
    novel_signatures = observed_signatures - expected_stop_signatures
    missing_expected_signatures = expected_stop_signatures - observed_signatures

    return {
        "total_stop_codon_count": stop_info["total_stop_codon_count"],
        "novel_stop_codon_count": len(novel_signatures),
        "missing_expected_stop_count": len(missing_expected_signatures),
        "stop_list": stop_list,
        "observed_signatures": observed_signatures,
        "novel_signatures": novel_signatures,
        "missing_expected_signatures": missing_expected_signatures,
    }


def simulate_one(sequence: str, mutation_rate: float, mutation_probabilities: dict,
                 orfs, expected_stop_signatures):
    """
    Simulate one descendant sequence from one parent.
    """
    seq = list(sequence)

    for i in range(len(seq)):
        if random.random() < mutation_rate:
            seq[i] = mutate_nucleotide(seq[i], mutation_probabilities)

    final_seq = "".join(seq)

    stop_summary = summarize_stops_against_expected(
        final_seq,
        orfs,
        expected_stop_signatures
    )
    stop_list = stop_summary["stop_list"]

    if stop_list:
        first = stop_list[0]
        stop_codon = first["stop_codon"]
        codon_start_1based = first["codon_start_1based"]
        codon_index_global_1based = first["codon_index_global_1based"]
        first_stop_orf_name = first["first_stop_orf_name"]
        codon_index_within_orf_1based = first["codon_index_within_orf_1based"]
    else:
        stop_codon = ""
        codon_start_1based = ""
        codon_index_global_1based = ""
        first_stop_orf_name = ""
        codon_index_within_orf_1based = ""

    return {
        "sequence": final_seq,
        "total_stop_codon_count": stop_summary["total_stop_codon_count"],
        "novel_stop_codon_count": stop_summary["novel_stop_codon_count"],
        "missing_expected_stop_count": stop_summary["missing_expected_stop_count"],
        "is_viable": (
            stop_summary["missing_expected_stop_count"] == 0
            and stop_summary["novel_stop_codon_count"] == 0
        ),
        "stop_codon": stop_codon,
        "codon_start_1based": codon_start_1based,
        "codon_index_global_1based": codon_index_global_1based,
        "first_stop_orf_name": first_stop_orf_name,
        "codon_index_within_orf_1based": codon_index_within_orf_1based,
        "mutated_nt_index_1based": "",
        "stop_signatures": stop_list,
    }


def simulate_parent_batch(args):
    """
    Simulate many descendants from one parent.
    """
    (
        parent_id,
        parent_seq,
        n_sims,
        mutation_rate,
        mutation_probabilities,
        orfs,
        expected_stop_signatures,
        seed,
    ) = args

    random.seed(seed)
    variants = []

    for sim_idx in range(1, n_sims + 1):
        result = simulate_one(
            parent_seq,
            mutation_rate,
            mutation_probabilities,
            orfs,
            expected_stop_signatures,
        )
        result["simulation_index"] = sim_idx
        result["source_parent_id"] = parent_id
        variants.append(result)

    return variants
