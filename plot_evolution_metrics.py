#!/usr/bin/env python3

import argparse
import csv
import os
import matplotlib.pyplot as plt


PLOT_CHOICES = [
    "similarity",
    "dnds",
    "novel_stop_count",
    "structurally_valid",
    "novel_stop_proportion",
    "stop_breakdown",
    "viable_parent_count",
    "viable_parent_proportion",
]


def load_summary_tsv(path):
    """
    Load one generation summary TSV.
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


def set_optional_ymax(ymax):
    if ymax is not None:
        bottom, _ = plt.ylim()
        plt.ylim(bottom=bottom, top=ymax)


def plot_similarity(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]

        generations = [r["generation"] for r in rows]
        sim_ref1 = [r["mean_similarity_index_ref1"] for r in rows]
        sim_ref2 = [r["mean_similarity_index_ref2"] for r in rows]
        std_ref1 = [r["std_similarity_index_ref1"] for r in rows]
        std_ref2 = [r["std_similarity_index_ref2"] for r in rows]

        plt.errorbar(
            generations,
            sim_ref1,
            yerr=std_ref1,
            fmt='-o',
            capsize=3,
            label=f"{label} ref1",
        )
        plt.errorbar(
            generations,
            sim_ref2,
            yerr=std_ref2,
            fmt='--o',
            capsize=3,
            label=f"{label} ref2",
        )

    plt.xlabel("Generation")
    plt.ylabel("Similarity Index (cosine distance)")
    plt.title("Mean Similarity Index over Time")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_dnds(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]

        generations = [r["generation"] for r in rows]
        mean_dnds = [r["mean_dnds_parents"] for r in rows]
        std_dnds = [r["std_dnds_parents"] for r in rows]

        plt.errorbar(
            generations,
            mean_dnds,
            yerr=std_dnds,
            fmt='-o',
            capsize=3,
            label=label,
        )

    plt.xlabel("Generation")
    plt.ylabel("Mean dN/dS across sampled parents")
    plt.title("dN/dS over Time")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_novel_stop_counts(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        stop_counts = [r["novel_stop_variant_count"] for r in rows]
        plt.plot(generations, stop_counts, marker='o', label=label)

    plt.xlabel("Generation")
    plt.ylabel("Number of progeny with novel stop issues")
    plt.title("Novel Stop Variants over Time")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_structurally_valid_counts(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        valid_counts = [r["structurally_valid_variant_count"] for r in rows]
        plt.plot(generations, valid_counts, marker='o', label=label)

    plt.xlabel("Generation")
    plt.ylabel("Number of structurally valid progeny")
    plt.title("Structurally Valid Variants over Time")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_novel_stop_proportions(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        novel = [r["proportion_variants_with_novel_stops"] for r in rows]
        plt.plot(generations, novel, marker='o', label=label)

    plt.xlabel("Generation")
    plt.ylabel("Proportion of all generated progeny with novel stops")
    plt.title("Novel Stop Proportion over Time")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_stop_issue_breakdown(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        novel = [r["proportion_stop_variants_with_novel_stops"] for r in rows]
        original_only = [r["proportion_stop_variants_with_only_original_stops"] for r in rows]

        plt.plot(generations, novel, marker='o', label=f"{label} novel")
        plt.plot(generations, original_only, marker='x', linestyle='--', label=f"{label} original-only")

    plt.xlabel("Generation")
    plt.ylabel("Proportion of stop-issue progeny")
    plt.title("Breakdown of Stop-Issue Variants")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_viable_parent_count(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        viable_counts = [r["reproductive_parent_count"] for r in rows]
        plt.plot(generations, viable_counts, marker='o', label=label)

    plt.axhline(y=100, linestyle='--', linewidth=1, label="Maximum (100)")
    plt.xlabel("Generation")
    plt.ylabel("Viable parents (out of 100 sampled)")
    plt.title("Viable Parents per 100 Sampled")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_viable_parent_proportion(datasets, output_path, ymax=None):
    plt.figure(figsize=(8, 5))

    for ds in datasets:
        rows = ds["rows"]
        label = ds["label"]
        generations = [r["generation"] for r in rows]
        viable_props = [r["viable_parent_proportion"] for r in rows]
        plt.plot(generations, viable_props, marker='o', label=label)

    plt.axhline(y=1.0, linestyle='--', linewidth=1, label="Maximum (1.0)")
    plt.xlabel("Generation")
    plt.ylabel("Proportion viable (out of 100)")
    plt.title("Proportion of Viable Parents per Generation")
    set_optional_ymax(ymax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary-tsv",
        action="append",
        required=True,
        help="Path to a generation_summaries.tsv file. Repeat for multiple runs.",
    )
    parser.add_argument(
        "--label",
        action="append",
        required=True,
        help="Label for the corresponding summary TSV. Repeat in same order as --summary-tsv.",
    )
    parser.add_argument("--output-dir", required=True)

    parser.add_argument(
        "--plots",
        nargs="+",
        choices=PLOT_CHOICES,
        default=PLOT_CHOICES,
        help="Subset of plots to generate.",
    )

    parser.add_argument("--ymax-similarity", type=float, default=None)
    parser.add_argument("--ymax-dnds", type=float, default=None)
    parser.add_argument("--ymax-novel-stop-count", type=float, default=None)
    parser.add_argument("--ymax-structurally-valid", type=float, default=None)
    parser.add_argument("--ymax-novel-stop-proportion", type=float, default=None)
    parser.add_argument("--ymax-stop-breakdown", type=float, default=None)
    parser.add_argument("--ymax-viable-parent-count", type=float, default=None)
    parser.add_argument("--ymax-viable-parent-proportion", type=float, default=None)

    args = parser.parse_args()

    if len(args.summary_tsv) != len(args.label):
        raise ValueError("You must provide the same number of --summary-tsv and --label arguments.")

    os.makedirs(args.output_dir, exist_ok=True)

    datasets = []
    for path, label in zip(args.summary_tsv, args.label):
        datasets.append({
            "label": label,
            "rows": load_summary_tsv(path),
        })

    if "similarity" in args.plots:
        plot_similarity(
            datasets,
            os.path.join(args.output_dir, "mean_similarity_index_over_time.png"),
            ymax=args.ymax_similarity,
        )

    if "dnds" in args.plots:
        plot_dnds(
            datasets,
            os.path.join(args.output_dir, "dnds_over_time.png"),
            ymax=args.ymax_dnds,
        )

    if "novel_stop_count" in args.plots:
        plot_novel_stop_counts(
            datasets,
            os.path.join(args.output_dir, "novel_stop_variant_count_over_time.png"),
            ymax=args.ymax_novel_stop_count,
        )

    if "structurally_valid" in args.plots:
        plot_structurally_valid_counts(
            datasets,
            os.path.join(args.output_dir, "structurally_valid_variant_count_over_time.png"),
            ymax=args.ymax_structurally_valid,
        )

    if "novel_stop_proportion" in args.plots:
        plot_novel_stop_proportions(
            datasets,
            os.path.join(args.output_dir, "novel_stop_proportions_over_time.png"),
            ymax=args.ymax_novel_stop_proportion,
        )

    if "stop_breakdown" in args.plots:
        plot_stop_issue_breakdown(
            datasets,
            os.path.join(args.output_dir, "novel_stop_breakdown_over_time.png"),
            ymax=args.ymax_stop_breakdown,
        )

    if "viable_parent_count" in args.plots:
        plot_viable_parent_count(
            datasets,
            os.path.join(args.output_dir, "viable_parent_count_over_time.png"),
            ymax=args.ymax_viable_parent_count,
        )

    if "viable_parent_proportion" in args.plots:
        plot_viable_parent_proportion(
            datasets,
            os.path.join(args.output_dir, "viable_parent_proportion_over_time.png"),
            ymax=args.ymax_viable_parent_proportion,
        )

    print(f"Wrote selected plots to {args.output_dir}")


if __name__ == "__main__":
    main()
