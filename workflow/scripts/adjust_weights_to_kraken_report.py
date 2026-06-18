import argparse
import re
from collections import defaultdict
from pathlib import Path


def get_weight(seq_id: str) -> int:
    m = re.search(r"\((\d+)\)$", seq_id)
    if m:
        return int(m.group(1))
    return 1

def get_consensus_extra_hits(
    kraken_output: str | Path,
) -> list[tuple[int, int]]:
    hits = []

    with open(kraken_output, "r") as fh:
        for line in fh:
            if not line.strip():
                continue

            blocks = line.rstrip("\n").split("\t")

            status = blocks[0]
            seq_id = blocks[1]
            taxid = int(blocks[2])

            if status != "C":
                continue

            if not seq_id.startswith("consensus_"):
                continue

            weight = get_weight(seq_id)
            extra_weight = weight - 1

            if extra_weight <= 0:
                continue

            hits.append((taxid, extra_weight))

    return hits

def count_indent(name: str) -> int:
    return len(name) - len(name.lstrip(" "))

def parse_kraken_report(report_path: str | Path) -> list[dict]:
    nodes = []

    with open(report_path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue

            blocks = line.rstrip("\n").split("\t", 5)

            if len(blocks) != 6:
                raise ValueError(f"Bad report line: {line}")

            percent, clade_count, taxon_count, rank, taxid, name = blocks

            node = {
                "percent": float(percent),
                "clade_count": int(clade_count),
                "taxon_count": int(taxon_count),
                "rank": rank,
                "taxid": int(taxid),
                "name": name,
                "depth": count_indent(name),
                "parent": None,
                "children": [],
            }

            nodes.append(node)

    return nodes             

def build_tree(nodes: list[dict]) -> None:
    stack = []

    for i, node in enumerate(nodes):
        depth = node["depth"]

        while stack and stack[-1][1] >= depth:
            stack.pop()

        if stack:
            parent_i = stack[-1][0]
            node["parent"] = parent_i
            nodes[parent_i]["children"].append(i)

        stack.append((i, depth))

def add_extra_weights_to_report(
    nodes: list[dict],
    extra_hits: list[tuple[int, int]],
) -> None:
    taxid_to_index = {
        node["taxid"]: i
        for i, node in enumerate(nodes)
    }

    for taxid, extra_weight in extra_hits:
        if taxid not in taxid_to_index:
            print(f"WARNING: taxid {taxid} not found in report")
            continue

        node_i = taxid_to_index[taxid]

        nodes[node_i]["taxon_count"] += extra_weight

        while node_i is not None:
            nodes[node_i]["clade_count"] += extra_weight
            node_i = nodes[node_i]["parent"]

    total = sum(
        node["clade_count"]
        for node in nodes
        if node["parent"] is None
    )

    for node in nodes:
        node["percent"] = node["clade_count"] / total * 100 if total else 0.0

def write_kraken_report(nodes: list[dict], out_path: str | Path) -> None:
    with open(out_path, "w") as out:
        for node in nodes:
            out.write(
                f"{node['percent']:.2f}\t"
                f"{node['clade_count']}\t"
                f"{node['taxon_count']}\t"
                f"{node['rank']}\t"
                f"{node['taxid']}\t"
                f"{node['name']}\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken-output", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--out", required=True)
    
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    weighted_hits = get_consensus_extra_hits(args.kraken_output)

    nodes = parse_kraken_report(args.report)
    build_tree(nodes)
    add_extra_weights_to_report(nodes, weighted_hits)
    write_kraken_report(nodes, args.out)


if __name__ == "__main__":
    main()
