#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd


RANK_LEVELS = {
    "D": 0,
    "K": 1,
    "P": 2,
    "C": 3,
    "O": 4,
    "F": 5,
    "G": 6,
    "S": 7,
}


def read_kraken_report_with_lineage(path: Path) -> pd.DataFrame:
    rows = []
    lineage = {
        "D": None,
        "K": None,
        "P": None,
        "C": None,
        "O": None,
        "F": None,
        "G": None,
        "S": None,
    }

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 6:
                continue

            percentage = float(parts[0].strip())
            reads_clade = int(parts[1].strip())
            reads_direct = int(parts[2].strip())
            rank = parts[3].strip()
            taxid = parts[4].strip()

            raw_name = parts[5]
            name = raw_name.strip()

            if rank in RANK_LEVELS:
                current_level = RANK_LEVELS[rank]

                for r, level in RANK_LEVELS.items():
                    if level >= current_level:
                        lineage[r] = None

                lineage[rank] = name

            rows.append({
                "percentage": percentage,
                "reads_clade": reads_clade,
                "reads_direct": reads_direct,
                "rank": rank,
                "taxid": taxid,
                "name": name,
                "domain": lineage["D"],
                "kingdom": lineage["K"],
                "phylum": lineage["P"],
                "class": lineage["C"],
                "order": lineage["O"],
                "family": lineage["F"],
                "genus": lineage["G"],
                "species": lineage["S"],
            })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--kraken_dir",
        required=True,
        help="Directory containing sample/report.tsv files"
    )
    parser.add_argument(
        "--rank",
        default="F",
        help="Taxonomic rank to extract: D, K, P, C, O, F, G, S"
    )
    parser.add_argument(
        "--min_reads",
        type=int,
        default=0,
        help="Minimum clade reads to keep a taxon in a sample"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file"
    )

    args = parser.parse_args()

    kraken_dir = Path(args.kraken_dir)
    report_files = sorted(kraken_dir.glob("*/report.tsv"))

    if not report_files:
        raise FileNotFoundError(f"No report.tsv files found in {kraken_dir}")

    all_tables = []
    metadata_tables = []

    for report in report_files:
        sample = report.parent.name

        df = read_kraken_report_with_lineage(report)

        df = df[df["rank"] == args.rank].copy()
        df = df[df["reads_clade"] >= args.min_reads].copy()

        meta_cols = [
            "name", "taxid", "domain", "kingdom",
            "phylum", "class", "order"
        ]

        metadata_tables.append(df[meta_cols].drop_duplicates("name"))

        sample_table = df[["name", "reads_clade"]].copy()
        sample_table = sample_table.rename(columns={"reads_clade": sample})
        sample_table = sample_table.groupby("name", as_index=False)[sample].sum()

        all_tables.append(sample_table)

    matrix = all_tables[0]

    for table in all_tables[1:]:
        matrix = matrix.merge(table, on="name", how="outer")

    matrix = matrix.fillna(0)

    sample_cols = [col for col in matrix.columns if col != "name"]
    matrix[sample_cols] = matrix[sample_cols].astype(int)

    metadata = pd.concat(metadata_tables, ignore_index=True)
    metadata = metadata.drop_duplicates("name")

    matrix = metadata.merge(matrix, on="name", how="right")

    matrix["total"] = matrix[sample_cols].sum(axis=1)
    matrix["present_in_samples"] = (matrix[sample_cols] > 0).sum(axis=1)

    front_cols = [
        "name", "taxid", "domain", "kingdom",
        "phylum", "class", "order",
        "present_in_samples", "total"
    ]

    matrix = matrix[front_cols + sample_cols]
    matrix = matrix.sort_values("total", ascending=False)

    matrix.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()