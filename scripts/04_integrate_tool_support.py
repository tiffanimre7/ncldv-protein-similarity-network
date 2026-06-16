"""
Integrate BLASTP, MMseqs2, and SWIPE support into consensus edge tables.

Usage:
    python scripts/04_integrate_tool_support.py

Required input files:
    - data/edges/<collapse_rule>/<threshold>/blastp.edges.tsv
    - data/edges/<collapse_rule>/<threshold>/mmseqs.edges.tsv
    - data/edges/<collapse_rule>/<threshold>/swipe.edges.tsv

These files are produced by scripts/02_build_edge_tables.py. The script combines the
three per-tool edge tables for each collapse rule and threshold into one integrated
network and writes a QC summary TSV.

This script reads the per-tool edge tables from scripts/02_build_edge_tables.py,
merges them into integrated edge tables, annotates tool support and support patterns,
and writes both the integrated tables and a QC summary TSV.
"""
from __future__ import annotations

from pathlib import Path
from functools import reduce
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

EDGE_ROOT = PROJECT_ROOT / "data" / "edges"
QC_OUTPUT_DIR = PROJECT_ROOT / "data" / "quality_checks"

TOOLS = ["blastp", "mmseqs", "swipe"]

COLLAPSE_RULES = [
    "best_bitscore",
    "best_evalue",
]

THRESHOLDS = [
    "raw",
    "evalue_1e-3",
    "evalue_1e-5",
    "evalue_1e-10",
]


def load_tool_edges(edge_file: Path, tool: str) -> pd.DataFrame:
    """Load one per-tool edge table and rename score columns so they remain traceable after cross-tool integration"""
    df = pd.read_csv(edge_file, sep="\t")

    required_cols = {
        "source",
        "target",
        "pident",
        "length",
        "evalue",
        "bitscore",
        "n_supporting_hits",
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"{edge_file} is missing required columns: {missing_cols}")

    df = df[
        [
            "source",
            "target",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "n_supporting_hits",
        ]
    ].copy()

    df = df.rename(
        columns={
            "pident": f"{tool}_pident",
            "length": f"{tool}_length",
            "evalue": f"{tool}_evalue",
            "bitscore": f"{tool}_bitscore",
            "n_supporting_hits": f"{tool}_supporting_hits",
        }
    )

    df[f"{tool}_supported"] = True

    return df


def merge_tool_edges(edge_tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Outer-merge edge tables from different tools by source-target pair"""
    return reduce(
        lambda left, right: pd.merge(left, right, on=["source", "target"], how="outer"),
        edge_tables,
    )


def add_support_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add n_tools, tools, support_pattern and combined score descriptors"""
    df = df.copy()

    for tool in TOOLS:
        supported_col = f"{tool}_supported"
        hits_col = f"{tool}_supporting_hits"

        if supported_col not in df.columns:
            df[supported_col] = False

        df[supported_col] = df[supported_col].fillna(False).astype(bool)

        if hits_col in df.columns:
            df[hits_col] = df[hits_col].fillna(0).astype(int)

    supported_cols = [f"{tool}_supported" for tool in TOOLS]
    df["n_tools"] = df[supported_cols].sum(axis=1)

    def make_tools(row: pd.Series) -> str:
        return ";".join([tool for tool in TOOLS if row[f"{tool}_supported"]])

    df["tools"] = df.apply(make_tools, axis=1)

    df["support_pattern"] = df["tools"].str.replace(";", "_", regex=False)

    pident_cols = [f"{tool}_pident" for tool in TOOLS if f"{tool}_pident" in df.columns]
    length_cols = [f"{tool}_length" for tool in TOOLS if f"{tool}_length" in df.columns]
    evalue_cols = [f"{tool}_evalue" for tool in TOOLS if f"{tool}_evalue" in df.columns]
    bitscore_cols = [f"{tool}_bitscore" for tool in TOOLS if f"{tool}_bitscore" in df.columns]
    hits_cols = [f"{tool}_supporting_hits" for tool in TOOLS if f"{tool}_supporting_hits" in df.columns]

    df["pident_max"] = df[pident_cols].max(axis=1)
    df["pident_mean"] = df[pident_cols].mean(axis=1)

    df["align_len_max"] = df[length_cols].max(axis=1)
    df["align_len_mean"] = df[length_cols].mean(axis=1)

    df["evalue_min"] = df[evalue_cols].min(axis=1)
    df["bitscore_max"] = df[bitscore_cols].max(axis=1)

    df["total_supporting_hits"] = df[hits_cols].sum(axis=1)

    front_cols = [
        "source",
        "target",
        "n_tools",
        "tools",
        "support_pattern",
        "blastp_supported",
        "mmseqs_supported",
        "swipe_supported",
        "pident_max",
        "pident_mean",
        "align_len_max",
        "align_len_mean",
        "evalue_min",
        "bitscore_max",
        "total_supporting_hits",
    ]

    remaining_cols = [col for col in df.columns if col not in front_cols]

    return df[front_cols + remaining_cols]


def summarize_integrated_edges(
    df: pd.DataFrame,
    collapse_rule: str,
    threshold: str,
    output_file: Path,
) -> dict:
    """Create one QC row for one integrated network"""
    support_counts = df["n_tools"].value_counts().to_dict()
    pattern_counts = df["support_pattern"].value_counts().to_dict()

    return {
        "collapse_rule": collapse_rule,
        "threshold": threshold,
        "output_file": str(output_file.relative_to(PROJECT_ROOT)),
        "integrated_edges": len(df),
        "edges_supported_by_1_tool": support_counts.get(1, 0),
        "edges_supported_by_2_tools": support_counts.get(2, 0),
        "edges_supported_by_3_tools": support_counts.get(3, 0),
        "blastp_edges": int(df["blastp_supported"].sum()),
        "mmseqs_edges": int(df["mmseqs_supported"].sum()),
        "swipe_edges": int(df["swipe_supported"].sum()),
        "blastp_only": pattern_counts.get("blastp", 0),
        "mmseqs_only": pattern_counts.get("mmseqs", 0),
        "swipe_only": pattern_counts.get("swipe", 0),
        "blastp_mmseqs": pattern_counts.get("blastp_mmseqs", 0),
        "blastp_swipe": pattern_counts.get("blastp_swipe", 0),
        "mmseqs_swipe": pattern_counts.get("mmseqs_swipe", 0),
        "blastp_mmseqs_swipe": pattern_counts.get("blastp_mmseqs_swipe", 0),
        "median_pident_max": df["pident_max"].median(),
        "median_align_len_max": df["align_len_max"].median(),
        "median_bitscore_max": df["bitscore_max"].median(),
        "min_evalue": df["evalue_min"].min(),
        "max_evalue": df["evalue_min"].max(),
    }


def process_combination(collapse_rule: str, threshold: str) -> dict | None:
    """Integrate BLASTP, MMseqs2 and SWIPE edge tables for one collapse rule and one threshold"""
    print(f"\nIntegrating {collapse_rule} | {threshold}")

    edge_tables = []

    for tool in TOOLS:
        edge_file = EDGE_ROOT / collapse_rule / threshold / f"{tool}.edges.tsv"

        if not edge_file.exists():
            print(f"  Missing file, skipping: {edge_file}")
            return None

        print(f"  Reading {tool}: {edge_file}")
        edge_tables.append(load_tool_edges(edge_file, tool))

    integrated = merge_tool_edges(edge_tables)
    integrated = add_support_columns(integrated)

    output_dir = EDGE_ROOT / "integrated" / collapse_rule / threshold
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "integrated_edges.tsv"
    integrated.to_csv(output_file, sep="\t", index=False)

    print(f"  Integrated edges written: {len(integrated):,}")
    print(f"  Output: {output_file}")

    pattern_output = output_dir / "support_pattern_counts.tsv"
    (
        integrated["support_pattern"]
        .value_counts()
        .rename_axis("support_pattern")
        .reset_index(name="edge_count")
        .to_csv(pattern_output, sep="\t", index=False)
    )

    return summarize_integrated_edges(
        df=integrated,
        collapse_rule=collapse_rule,
        threshold=threshold,
        output_file=output_file,
    )


def main() -> None:
    QC_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    qc_rows = []

    for collapse_rule in COLLAPSE_RULES:
        for threshold in THRESHOLDS:
            qc_row = process_combination(collapse_rule, threshold)

            if qc_row is not None:
                qc_rows.append(qc_row)

    qc_df = pd.DataFrame(qc_rows)

    qc_output = QC_OUTPUT_DIR / "integrated_tool_support_qc.tsv"
    qc_df.to_csv(qc_output, sep="\t", index=False)

    print("\nIntegrated tool support QC:")
    print(qc_df.to_string(index=False))
    print(f"\nQC table written to: {qc_output}")


if __name__ == "__main__":
    main()