"""
Build edge tables from standardized alignments and write a QC summary.

Usage:
    python scripts/02_build_edge_tables.py

Required input files:
    - data/processed/blastp.standardized.tsv
    - data/processed/mmseqs.standardized.tsv
    - data/processed/swipe.standardized.tsv

These files are produced by scripts/01_standardize_alignments.py. If you are using
the scripts independently, make sure those standardized TSV files exist before
running this step.

This script reads the standardized alignment files for BLASTP, MMseqs2, and SWIPE,
converts query-subject pairs into undirected protein-pair edges, applies multiple
E-value thresholds, collapses duplicate hits by the chosen rule, and writes edge
tables plus a summary TSV under data/quality_checks/.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUTS = {
    "blastp": PROJECT_ROOT / "data" / "processed" / "blastp.standardized.tsv",
    "mmseqs": PROJECT_ROOT / "data" / "processed" / "mmseqs.standardized.tsv",
    "swipe": PROJECT_ROOT / "data" / "processed" / "swipe.standardized.tsv",
}

EDGE_OUTPUT_DIR = PROJECT_ROOT / "data" / "edges"
QC_OUTPUT_DIR = PROJECT_ROOT / "data" / "quality_checks"

THRESHOLDS = {
    "raw": None,
    "evalue_1e-3": 1e-3,
    "evalue_1e-5": 1e-5,
    "evalue_1e-10": 1e-10,
} # the more the better, the more stringent the threshold, the less edges will be kept

COLLAPSE_RULES = [
    "best_bitscore", # we take the bitscore parameter to collapse the edges, the higher the better
    "best_evalue",
]


def add_undirected_edge_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Convert directed query-subject pairs into undirected protein-protein edges"""
    df = df.copy()

    df["source"] = np.where(
        df["qseqid"] <= df["sseqid"],
        df["qseqid"],
        df["sseqid"],
    )

    df["target"] = np.where(
        df["qseqid"] <= df["sseqid"],
        df["sseqid"],
        df["qseqid"],
    )

    return df


def apply_threshold(df: pd.DataFrame, threshold_value: float | None) -> pd.DataFrame:
    """Apply an E-value threshold to filter alignment hits, if threshold_value is None, keep the raw reported hits"""
    if threshold_value is None:
        return df.copy()

    return df[df["evalue"] <= threshold_value].copy()


def collapse_edges(df: pd.DataFrame, collapse_rule: str) -> pd.DataFrame:
    """
    Collapse multiple alignment hits between the same protein pair into one edge

    best_bitscore:
        keep the hit with the highest bit score

    best_evalue:
        keep the hit with the lowest E-value
    """
    group_cols = ["source", "target"]

    supporting_hits = (
        df.groupby(group_cols)
        .size()
        .reset_index(name="n_supporting_hits")
    )

    if collapse_rule == "best_bitscore":
        best_idx = df.groupby(group_cols)["bitscore"].idxmax()
    elif collapse_rule == "best_evalue":
        best_idx = df.groupby(group_cols)["evalue"].idxmin()
    else:
        raise ValueError(f"Unknown collapse rule: {collapse_rule}")

    edges = df.loc[best_idx].copy()

    edges = edges.merge(
        supporting_hits,
        on=group_cols,
        how="left",
    )

    edge_cols = [
        "source",
        "target",
        "tool",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "evalue",
        "bitscore",
        "n_supporting_hits",
    ]

    return edges[edge_cols].sort_values(["source", "target"])


def summarize_edges(
    tool: str,
    threshold_name: str,
    threshold_value: float | None,
    collapse_rule: str,
    rows_before_threshold: int,
    rows_after_threshold: int,
    edges: pd.DataFrame,
) -> dict:
    """Create one QC summary row for one tool / threshold / collapse-rule combination"""
    return {
        "tool": tool,
        "threshold": threshold_name,
        "threshold_value": threshold_value if threshold_value is not None else "raw",
        "collapse_rule": collapse_rule,
        "alignment_rows_before_threshold": rows_before_threshold,
        "alignment_rows_after_threshold": rows_after_threshold,
        "unique_edges": len(edges),
        "collapsed_duplicate_hits": rows_after_threshold - len(edges),
        "min_pident": edges["pident"].min(),
        "median_pident": edges["pident"].median(),
        "max_pident": edges["pident"].max(),
        "min_length": edges["length"].min(),
        "median_length": edges["length"].median(),
        "max_length": edges["length"].max(),
        "min_evalue": edges["evalue"].min(),
        "max_evalue": edges["evalue"].max(),
        "min_bitscore": edges["bitscore"].min(),
        "median_bitscore": edges["bitscore"].median(),
        "max_bitscore": edges["bitscore"].max(),
        "min_supporting_hits": edges["n_supporting_hits"].min(),
        "max_supporting_hits": edges["n_supporting_hits"].max(),
    }


def process_tool(tool: str, input_path: Path) -> list[dict]:
    """Process one standardized alignment file and create edge tables"""
    print(f"\nReading standardized alignments for {tool}:") # so we can see which tool is being processed real time
    print(input_path)

    df = pd.read_csv(input_path, sep="\t")

    required_cols = {
        "tool", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"{tool} is missing required columns: {missing_cols}")

    df = add_undirected_edge_columns(df)

    rows_before_threshold = len(df)

    qc_rows = []

    for threshold_name, threshold_value in THRESHOLDS.items():
        thresholded = apply_threshold(df, threshold_value)
        rows_after_threshold = len(thresholded)

        print(
            f"{tool} | {threshold_name}: "
            f"{rows_after_threshold:,} alignment rows after threshold"
        )

        if rows_after_threshold == 0:
            print(f"Skipping {tool} | {threshold_name}, no rows left.")
            continue

        for collapse_rule in COLLAPSE_RULES:
            edges = collapse_edges(thresholded, collapse_rule)

            output_dir = EDGE_OUTPUT_DIR / collapse_rule / threshold_name
            output_dir.mkdir(parents=True, exist_ok=True)

            output_path = output_dir / f"{tool}.edges.tsv"
            edges.to_csv(output_path, sep="\t", index=False)

            print(
                f"  {collapse_rule}: "
                f"{len(edges):,} unique edges written to {output_path}"
            )

            qc_rows.append(
                summarize_edges(
                    tool=tool,
                    threshold_name=threshold_name,
                    threshold_value=threshold_value,
                    collapse_rule=collapse_rule,
                    rows_before_threshold=rows_before_threshold,
                    rows_after_threshold=rows_after_threshold,
                    edges=edges,
                )
            )

    return qc_rows


def main() -> None:
    EDGE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_qc_rows = []

    for tool, input_path in INPUTS.items():
        all_qc_rows.extend(process_tool(tool, input_path))

    qc_df = pd.DataFrame(all_qc_rows)

    qc_output = QC_OUTPUT_DIR / "edge_table_construction_qc.tsv"
    qc_df.to_csv(qc_output, sep="\t", index=False)

    print("\nEdge table construction QC:")
    print(qc_df.to_string(index=False))
    print(f"\nQC table written to: {qc_output}")


if __name__ == "__main__":
    main()