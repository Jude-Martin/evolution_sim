#!/usr/bin/env python3

import argparse
import csv
import os
from collections import Counter, defaultdict

import matplotlib.pyplot as plt


def parse_label(path):
    base = os.path.basename(path)
    return os.path.splitext(base)[0]


def load_shaded_regions(path):
    """
    TSV columns:
      region_name    start_1based    end_1based
    """
    regions = []
    if path is None:
        return regions

    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_name", "start_1based", "end_1based"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(f"Shade-region TSV must contain columns: {sorted(required)}")

        for row in reader:
            start = int(row["start_1based"])
            end = int(row["end_1based"])
            if start < 1 or end < start:
                raise ValueError(f"Invalid shaded region: {row}")

            regions.append({
                "region_name": row["region_name"],
                "start_1based": start,
                "end_1based": end,
            })

    return regions


def load_novel_stop_events(paths, labels=None):
    """
    Load novel stop events from one or more gen_*_novel_stop_events.tsv files.
    """
    all_rows = []

    if labels is None:
        labels = [parse_label(p) for p in paths]

    if len(paths) != len(labels):
        raise ValueError("paths and labels must have same length")

    for path, label in zip(paths, labels):
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "codon_start_1based",
                "stop_codon",
                "generation",
                "simulation_index",
                "source_parent_id",
                "orf_name",
                "codon_index_within_orf_1based",
            }
            if not required.issubset(reader.fieldnames or set()):
                raise ValueError(f"{path} is missing required columns: {sorted(required)}")

            for row in reader:
                pos_int = int(row["codon_start_1based"])
                all_rows.append({
                    "source_label": label,
                    "generation": int(row["generation"]),
                    "simulation_index": int(row["simulation_index"]),
                    "source_parent_id": row["source_parent_id"],
                    "position_1based": pos_int,
                    "stop_codon": row["stop_codon"],
                    "orf_name": row["orf_name"],
                    "codon_index_within_orf_1based": int(row["codon_index_within_orf_1based"]),
                })

    return all_rows


def compute_hotspot_windows(stop_rows, sequence_length, window_width):
    half = window_width // 2
    positions = [r["position_1based"] for r in stop_rows]

    window_rows = []
    for center in range(1, sequence_length + 1):
        start = max(1, center - half)
        end = min(sequence_length, start + window_width - 1)
        start = max(1, end - window_width + 1)

        count = sum(1 for p in positions if start <= p <= end)
        window_rows.append({
            "window_start_1based": start,
            "window_end_1based": end,
            "window_center_1based": center,
            "stop_count": count,
        })

    return window_rows


def write_hotspot_tables(stop_rows, sequence_length, window_width, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    pos_counter = Counter(r["position_1based"] for r in stop_rows)

    position_counts_path = os.path.join(output_dir, "position_counts.tsv")
    with open(position_counts_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["position_1based", "stop_count"])
        for pos in range(1, sequence_length + 1):
            writer.writerow([pos, pos_counter.get(pos, 0)])

    codon_counts_path = os.path.join(output_dir, "position_codon_counts.tsv")
    codon_counter = defaultdict(Counter)
    for row in stop_rows:
        codon_counter[row["position_1based"]][row["stop_codon"]] += 1

    with open(codon_counts_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["position_1based", "TAA_count", "TAG_count", "TGA_count", "total_stop_count"])
        for pos in range(1, sequence_length + 1):
            taa = codon_counter[pos].get("TAA", 0)
            tag = codon_counter[pos].get("TAG", 0)
            tga = codon_counter[pos].get("TGA", 0)
            writer.writerow([pos, taa, tag, tga, taa + tag + tga])

    windows = compute_hotspot_windows(stop_rows, sequence_length, window_width)

    windows_path = os.path.join(output_dir, "hotspot_windows.tsv")
    with open(windows_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "window_start_1based",
                "window_end_1based",
                "window_center_1based",
                "stop_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(windows)

    sorted_path = os.path.join(output_dir, "hotspot_windows_sorted.tsv")
    with open(sorted_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "window_start_1based",
                "window_end_1based",
                "window_center_1based",
                "stop_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in sorted(windows, key=lambda x: x["stop_count"], reverse=True):
            writer.writerow(row)

    return windows


def add_shaded_regions_x(regions):
    for region in regions:
        plt.axvspan(region["start_1based"], region["end_1based"], alpha=0.15)


def add_shaded_regions_y(regions):
    for region in regions:
        plt.axhspan(region["start_1based"], region["end_1based"], alpha=0.15)


def plot_dotplot(stop_rows, sequence_length, shaded_regions, output_path):
    labels = sorted({r["source_label"] for r in stop_rows})
    label_to_y = {lab: i for i, lab in enumerate(labels)}

    plt.figure(figsize=(12, 6))
    add_shaded_regions_x(shaded_regions)

    xs = [r["position_1based"] for r in stop_rows]
    ys = [label_to_y[r["source_label"]] for r in stop_rows]

    plt.scatter(xs, ys, s=10)
    plt.xlim(1, sequence_length)
    plt.yticks(range(len(labels)), labels)
    plt.xlabel("Sequence position (1-based)")
    plt.ylabel("Input file / group")
    plt.title("Novel stop-codon occurrence dot plot")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_violin(stop_rows, sequence_length, shaded_regions, output_path):
    grouped = defaultdict(list)
    for row in stop_rows:
        grouped[row["source_label"]].append(row["position_1based"])

    labels = list(sorted(grouped.keys()))
    data = [grouped[label] for label in labels]

    plt.figure(figsize=(10, 6))
    add_shaded_regions_y(shaded_regions)

    if data:
        plt.violinplot(data, showmeans=False, showmedians=True, showextrema=True)

    plt.xticks(range(1, len(labels) + 1), labels, rotation=30, ha="right")
    plt.ylim(1, sequence_length)
    plt.ylabel("Sequence position (1-based)")
    plt.xlabel("Input file / group")
    plt.title("Novel stop-codon position distribution")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_moving_window(windows, shaded_regions, output_path):
    centers = [w["window_center_1based"] for w in windows]
    counts = [w["stop_count"] for w in windows]

    plt.figure(figsize=(12, 6))
    add_shaded_regions_x(shaded_regions)
    plt.plot(centers, counts)
    plt.xlabel("Sequence position (window center, 1-based)")
    plt.ylabel("Novel stop occurrences in window")
    plt.title("Moving-window novel stop hotspot profile")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Tabulate and plot novel stop-codon hotspots from gen_*_novel_stop_events.tsv files."
    )
    parser.add_argument(
        "--input-tsv",
        action="append",
        required=True,
        help="Input novel stop event TSV file(s). Repeat for multiple files."
    )
    parser.add_argument(
        "--label",
        action="append",
        default=None,
        help="Optional label(s) matching --input-tsv order."
    )
    parser.add_argument(
        "--sequence-length",
        type=int,
        required=True,
        help="Full sequence length."
    )
    parser.add_argument(
        "--window-width",
        type=int,
        default=50,
        help="Moving-window width. Default: 50"
    )
    parser.add_argument(
        "--shade-regions",
        default=None,
        help="Optional TSV with region_name/start_1based/end_1based for background shading."
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for plots and tables."
    )
    args = parser.parse_args()

    if args.label is not None and len(args.label) != len(args.input_tsv):
        raise ValueError("If provided, --label must be repeated the same number of times as --input-tsv.")

    os.makedirs(args.output_dir, exist_ok=True)

    shaded_regions = load_shaded_regions(args.shade_regions)
    stop_rows = load_novel_stop_events(args.input_tsv, args.label)

    if not stop_rows:
        raise ValueError("No novel stop events found in the supplied TSV file(s).")

    windows = write_hotspot_tables(
        stop_rows=stop_rows,
        sequence_length=args.sequence_length,
        window_width=args.window_width,
        output_dir=args.output_dir
    )

    plot_dotplot(
        stop_rows=stop_rows,
        sequence_length=args.sequence_length,
        shaded_regions=shaded_regions,
        output_path=os.path.join(args.output_dir, "stop_dotplot.png")
    )

    plot_violin(
        stop_rows=stop_rows,
        sequence_length=args.sequence_length,
        shaded_regions=shaded_regions,
        output_path=os.path.join(args.output_dir, "stop_violin.png")
    )

    plot_moving_window(
        windows=windows,
        shaded_regions=shaded_regions,
        output_path=os.path.join(args.output_dir, "stop_moving_window.png")
    )

    print(f"Wrote hotspot tables and plots to {args.output_dir}")


if __name__ == "__main__":
    main()
