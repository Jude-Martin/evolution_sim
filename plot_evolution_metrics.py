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
                "stop_variant_count": int(row["stop_variant_count"]),
                "reference_used": row["reference_used"],
            })
    return rows


def plot_similarity(rows, output_path):
    """
    Plot mean Similarity Index to both references over time.

    Note:
    Similarity Index here is cosine distance, so lower = more similar.
    """
    generations = [r["generation"] for r in rows]
    sim_ref1 = [r["mean_similarity_index_ref1"] for r in rows]
    sim_ref2 = [r["mean_similarity_index_ref2"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, sim_ref1, label="Mean Similarity Index to ref1")
    plt.plot(generations, sim_ref2, label="Mean Similarity Index to ref2")
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
    plt.plot(generations, stop_counts)
    plt.xlabel("Generation")
    plt.ylabel("Number of Stop-Containing Variants")
    plt.title("Stop Variants Generated over Time")
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

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
