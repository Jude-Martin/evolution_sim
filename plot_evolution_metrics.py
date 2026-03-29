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
            sampled_parent_count = int(row["sampled_parent_count"])
            reproductive_parent_count = int(row["reproductive_parent_count"])

            viable_parent_proportion = (
                reproductive_parent_count / sampled_parent_count
                if sampled_parent_count > 0 else 0.0
            )

            rows.append({
                "generation": int(row["generation"]),
                "sampled_parent_count": sampled_parent_count,
                "reproductive_parent_count": reproductive_parent_count,
                "viable_parent_proportion": viable_parent_proportion,
                "mean_similarity_index_ref1": float(row["mean_similarity_index_ref1"]),
                "mean_similarity_index_ref2": float(row["mean_similarity_index_ref2"]),
                "std_similarity_index_ref1": float(row["std_similarity_index_ref1"]),
                "std_similarity_index_ref2": float(row["std_similarity_index_ref2"]),
                "mean_dnds_parents": float(row["mean_dnds_parents"]),
                "std_dnds_parents": float(row["std_dnds_parents"]),
                "novel_stop_variant_count": int(row["novel_stop_variant_count"]),
                "structurally_valid_variant_count": int(row["structurally_valid_variant_count"]),
                "proportion_variants_with_novel_stops": float(row["proportion_variants_with_novel_stops"]),
                "proportion_variants_with_only_original_stops": float(row["proportion_variants_with_only_original_stops"]),
                "proportion_stop_variants_with_novel_stops": float(row["proportion_stop_variants_with_novel_stops"]),
                "proportion_stop_variants_with_only_original_stops": float(row["proportion_stop_variants_with_only_original_stops"]),
                "reference_used": row["reference_used"],
            })
    return rows


def plot_similarity(rows, output_path):
    """
    Plot mean Similarity Index to both references over time with SD error bars.
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


def plot_dnds(rows, output_path):
    """
    Plot mean dN/dS across sampled parents with SD error bars.
    """
    generations = [r["generation"] for r in rows]
    mean_dnds = [r["mean_dnds_parents"] for r in rows]
    std_dnds = [r["std_dnds_parents"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.errorbar(generations, mean_dnds, yerr=std_dnds, fmt='-o', capsize=3)
    plt.xlabel("Generation")
    plt.ylabel("Mean dN/dS across sampled parents")
    plt.title("dN/dS over Time")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_novel_stop_counts(rows, output_path):
    """
    Plot the number of progeny per generation that contain novel stop-position issues.
    """
    generations = [r["generation"] for r in rows]
    stop_counts = [r["novel_stop_variant_count"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, stop_counts, marker='o')
    plt.xlabel("Generation")
    plt.ylabel("Number of progeny with novel stop issues")
    plt.title("Novel Stop Variants over Time")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_structurally_valid_counts(rows, output_path):
    """
    Plot the number of structurally valid progeny per generation.
    """
    generations = [r["generation"] for r in rows]
    valid_counts = [r["structurally_valid_variant_count"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, valid_counts, marker='o')
    plt.xlabel("Generation")
    plt.ylabel("Number of structurally valid progeny")
    plt.title("Structurally Valid Variants over Time")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_novel_stop_proportions(rows, output_path):
    """
    Plot the proportion of all generated progeny that contain:
    - novel stop positions
    - only original stop-position issues
    """
    generations = [r["generation"] for r in rows]
    novel = [r["proportion_variants_with_novel_stops"] for r in rows]
    original_only = [r["proportion_variants_with_only_original_stops"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, novel, marker='o', label="Novel stop issues")
    plt.plot(generations, original_only, marker='o', label="Only original stop issues")
    plt.xlabel("Generation")
    plt.ylabel("Proportion of all generated progeny")
    plt.title("Novel Stop Composition over Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_stop_issue_breakdown(rows, output_path):
    """
    Plot the proportion among stop-issue progeny only:
    - novel stop issues
    - only original stop issues
    """
    generations = [r["generation"] for r in rows]
    novel = [r["proportion_stop_variants_with_novel_stops"] for r in rows]
    original_only = [r["proportion_stop_variants_with_only_original_stops"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, novel, marker='o', label="Novel stop issues")
    plt.plot(generations, original_only, marker='o', label="Only original stop issues")
    plt.xlabel("Generation")
    plt.ylabel("Proportion of stop-issue progeny")
    plt.title("Breakdown of Stop-Issue Variants")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_viable_parent_count(rows, output_path):
    """
    Plot the number of viable parents out of 100 sampled per generation.
    Includes a reference line at 100.
    """
    generations = [r["generation"] for r in rows]
    viable_counts = [r["reproductive_parent_count"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, viable_counts, marker='o', label="Viable parents")
    plt.axhline(y=100, linestyle='--', linewidth=1, label="Maximum (100)")
    plt.xlabel("Generation")
    plt.ylabel("Viable parents (out of 100 sampled)")
    plt.title("Viable Parents per 100 Sampled")
    plt.ylim(0, 105)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_viable_parent_proportion(rows, output_path):
    """
    Plot the proportion of viable parents per generation.
    Includes a reference line at 1.0.
    """
    generations = [r["generation"] for r in rows]
    viable_props = [r["viable_parent_proportion"] for r in rows]

    plt.figure(figsize=(8, 5))
    plt.plot(generations, viable_props, marker='o', label="Viable proportion")
    plt.axhline(y=1.0, linestyle='--', linewidth=1, label="Maximum (1.0)")
    plt.xlabel("Generation")
    plt.ylabel("Proportion viable (out of 100)")
    plt.title("Proportion of Viable Parents per Generation")
    plt.ylim(0, 1.05)
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

    plot_similarity(rows, os.path.join(args.output_dir, "mean_similarity_index_over_time.png"))
    plot_dnds(rows, os.path.join(args.output_dir, "dnds_over_time.png"))
    plot_novel_stop_counts(rows, os.path.join(args.output_dir, "novel_stop_variant_count_over_time.png"))
    plot_structurally_valid_counts(rows, os.path.join(args.output_dir, "structurally_valid_variant_count_over_time.png"))
    plot_novel_stop_proportions(rows, os.path.join(args.output_dir, "novel_stop_proportions_over_time.png"))
    plot_stop_issue_breakdown(rows, os.path.join(args.output_dir, "novel_stop_breakdown_over_time.png"))
    plot_viable_parent_count(rows, os.path.join(args.output_dir, "viable_parent_count_over_time.png"))
    plot_viable_parent_proportion(rows, os.path.join(args.output_dir, "viable_parent_proportion_over_time.png"))

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
