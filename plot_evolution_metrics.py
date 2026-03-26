import argparse
import csv
import os
import matplotlib.pyplot as plt


def load_summary_tsv(path):
    """
    Load generation summary metrics from TSV.
    """
    rows = []
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({
                "generation": int(row["generation"]),
                "mean_similarity_index_ref1": float(row["mean_similarity_index_ref1"]),
                "mean_similarity_index_ref2": float(row["mean_similarity_index_ref2"]),
                "std_similarity_index_ref1": float(row["std_similarity_index_ref1"]),
                "std_similarity_index_ref2": float(row["std_similarity_index_ref2"]),
                "stop_variant_count": int(row["stop_variant_count"]),
                "proportion_variants_with_novel_stops": float(row["proportion_variants_with_novel_stops"]),
                "proportion_variants_with_only_original_stops": float(row["proportion_variants_with_only_original_stops"]),
                "reference_used": row["reference_used"],
            })
    return rows


def plot_similarity(rows, output_path):
    """
    Plot mean Similarity Index to both references over time with
    standard-deviation error bars.

    Similarity Index here is cosine distance, so lower = more similar.
    """
    generations = [r["generation"] for r in rows]

    sim_ref1 = [r["mean_similarity_index_ref1"] for r in rows]
    sim_ref2 = [r["mean_similarity_index_ref2"] for r in rows]

    std_ref1 = [r["std_similarity_index_ref1"] for r in rows]
    std_ref2 = [r["std_similarity_index_ref2"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.errorbar(generations, sim_ref1, yerr=std_ref1, fmt='-o', capsize=3, label="Mean Similarity Index to ref1")
    plt.errorbar(generations, sim_ref2, yerr=std_ref2, fmt='-o', capsize=3, label="Mean Similarity Index to ref2")
    plt.xlabel("Generation")
    plt.ylabel("Similarity Index (cosine distance)")
    plt.title("Mean Similarity Index over Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_stop_counts(rows, output_path):
    """
    Plot number of stop-containing variants generated per generation.
    """
    generations = [r["generation"] for r in rows]
    stop_counts = [r["stop_variant_count"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, stop_counts, marker='o')
    plt.xlabel("Generation")
    plt.ylabel("Number of Stop-Containing Variants")
    plt.title("Stop Variants Generated over Time")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_stop_proportions(rows, output_path):
    """
    Plot proportion of each generation's variants that contain:
    - novel stop codons not present in the original sequence
    - only stop codons already present in the original sequence
    """
    generations = [r["generation"] for r in rows]
    novel = [r["proportion_variants_with_novel_stops"] for r in rows]
    original_only = [r["proportion_variants_with_only_original_stops"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, novel, marker='o', label="Novel stop codons")
    plt.plot(generations, original_only, marker='o', label="Only original stop codons")
    plt.xlabel("Generation")
    plt.ylabel("Proportion of all generated variants")
    plt.title("Stop-Codon Composition over Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    rows = load_summary_tsv(args.summary_tsv)

    plot_similarity(
        rows,
        os.path.join(args.output_dir, "mean_similarity_index_over_time.png")
    )
    plot_stop_counts(
        rows,
        os.path.join(args.output_dir, "stop_variant_count_over_time.png")
    )
    plot_stop_proportions(
        rows,
        os.path.join(args.output_dir, "stop_proportions_over_time.png")
    )

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
