import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot contaminant hit coverage along normalized read body."
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Length and hits joined table",
    )

    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output plot file",
    )

    parser.add_argument(
        "--bins", "-b",
        type=int,
        default=100,
        help="Number of bins along normalized read body",
    )

    return parser.parse_args()


def calculate_coverage(n_bins: int, df: pd.DataFrame) -> list[dict]:
    rows = []

    for pattern_name, sub in df.groupby("patternName"):
        coverage = np.zeros(n_bins, dtype=int)

        for _, row in sub.iterrows():
            start_bin = int(np.floor(row["rel_start"] / 100 * n_bins))
            end_bin = int(np.ceil(row["rel_end"] / 100 * n_bins))

            start_bin = max(0, min(n_bins - 1, start_bin))
            end_bin = max(start_bin + 1, min(n_bins, end_bin))

            coverage[start_bin:end_bin] += 1

        for i, cov in enumerate(coverage):
            rows.append({
                "patternName": pattern_name,
                "bin": i,
                "position_percent": (i + 0.5) / n_bins * 100,
                "coverage": cov,
            })

    return rows


def plot_empty(outfile: str, sample_name: str) -> None:
    plt.figure(figsize=(12, 5))
    plt.text(
        0.5,
        0.5,
        "No contaminant hits found",
        ha="center",
        va="center",
        transform=plt.gca().transAxes,
        fontsize=14,
    )
    plt.xlabel("Relative position in read (%)")
    plt.ylabel("Hit coverage")
    plt.title(f"{sample_name}: contaminant body coverage")
    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def plot_cov(coverage: list[dict], outfile: str, sample_name: str) -> None:
    cov_df = pd.DataFrame(coverage)

    if cov_df.empty:
        plot_empty(outfile, sample_name)
        return

    plt.figure(figsize=(12, 5))

    for pattern_name, sub in cov_df.groupby("patternName"):
        plt.plot(
            sub["position_percent"],
            sub["coverage"],
            label=pattern_name,
            linewidth=1.5,
        )

    plt.xlabel("Relative position in read (%)")
    plt.ylabel("Hit coverage")
    plt.title(f"{sample_name}: contaminant body coverage")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def main() -> None:
    args = parse_args()

    input_file = args.input
    out_file = args.output
    sample = Path(input_file).stem

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)

    try:
        df = pd.read_csv(input_file, sep="\t")
    except pd.errors.EmptyDataError:
        plot_empty(out_file, sample)
        return

    if df.empty:
        plot_empty(out_file, sample)
        return

    required_cols = {"patternName", "rel_start", "rel_end"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {input_file}: {missing}")

    df["rel_start"] = pd.to_numeric(df["rel_start"], errors="coerce")
    df["rel_end"] = pd.to_numeric(df["rel_end"], errors="coerce")
    df = df.dropna(subset=["rel_start", "rel_end"])

    if df.empty:
        plot_empty(out_file, sample)
        return

    coverage = calculate_coverage(args.bins, df)
    plot_cov(coverage, out_file, sample)


if __name__ == "__main__":
    main()