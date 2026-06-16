"""
Export consensus edge tables and Cytoscape-ready network files.

Usage:
    python scripts/06_export_consensus_and_cytoscape.py

Required input files:
    - data/edges/integrated/best_bitscore/evalue_1e-5/integrated_edges.tsv

This script only works on the primary integrated network at E-value <= 1e-5.
That integrated edge table is produced by scripts/04_integrate_tool_support.py.

This script reads the primary integrated edge table at E-value <= 1e-5, writes
consensus/all-edge subsets, creates Cytoscape-friendly edge and node tables for the
full consensus_3tools network and top-bit-score visualizations, and writes a QC
summary TSV.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INTEGRATED_INPUT = (
    PROJECT_ROOT
    / "data"
    / "edges"
    / "integrated"
    / "best_bitscore"
    / "evalue_1e-5"
    / "integrated_edges.tsv"
)

CONSENSUS_OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "edges"
    / "consensus"
    / "best_bitscore"
    / "evalue_1e-5"
)

CYTOSCAPE_OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "cytoscape"
    / "best_bitscore"
    / "evalue_1e-5"
)

QC_OUTPUT_DIR = PROJECT_ROOT / "data" / "quality_checks"

NETWORK_FILTERS = {
    "all_edges": "n_tools >= 1",
    "consensus_2plus": "n_tools >= 2",
    "consensus_3tools": "n_tools == 3",
}

TOP_N_VALUES = [10000, 50000]

KEEP_COLS = [
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


def load_integrated_edges() -> pd.DataFrame:
    """Load the primary integrated edge table"""
    print(f"Reading integrated edge table:\n{INTEGRATED_INPUT}")

    df = pd.read_csv(INTEGRATED_INPUT, sep="\t", usecols=KEEP_COLS)

    required_cols = {
        "source",
        "target",
        "n_tools",
        "tools",
        "support_pattern",
        "pident_max",
        "align_len_max",
        "evalue_min",
        "bitscore_max",
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return df


def filter_network(df: pd.DataFrame, filter_name: str) -> pd.DataFrame:
    """Create all_edges, consensus_2plus, or consensus_3tools subset"""
    if filter_name == "all_edges":
        return df[df["n_tools"] >= 1].copy()

    if filter_name == "consensus_2plus":
        return df[df["n_tools"] >= 2].copy()

    if filter_name == "consensus_3tools":
        return df[df["n_tools"] == 3].copy()

    raise ValueError(f"Unknown filter_name: {filter_name}")


def add_cytoscape_columns(edges: pd.DataFrame) -> pd.DataFrame:
    """Add Cytoscape-friendly columns"""
    cyt = edges.copy()

    cyt["interaction"] = "similarity"

    # Main edge weight for visualization
    cyt["edge_weight"] = cyt["bitscore_max"]

    # Avoid infinite values, when e-value is reported as 0.0
    clipped_evalue = cyt["evalue_min"].clip(lower=1e-300)
    cyt["neg_log10_evalue"] = -np.log10(clipped_evalue)

    cyt["support_class"] = np.where(
        cyt["n_tools"] == 3,
        "three_tools",
        np.where(cyt["n_tools"] == 2, "two_tools", "one_tool"),
    )

    cyt_cols = [
        "source",
        "target",
        "interaction",
        "edge_weight",
        "support_class",
        "n_tools",
        "tools",
        "support_pattern",
        "pident_max",
        "pident_mean",
        "align_len_max",
        "align_len_mean",
        "evalue_min",
        "neg_log10_evalue",
        "bitscore_max",
        "total_supporting_hits",
        "blastp_supported",
        "mmseqs_supported",
        "swipe_supported",
    ]

    return cyt[cyt_cols]


def make_node_table(cyt_edges: pd.DataFrame) -> pd.DataFrame:
    """Create a simple Cytoscape node attribute table"""
    all_nodes = pd.concat(
        [cyt_edges["source"], cyt_edges["target"]],
        ignore_index=True,
    )

    node_degree = (
        all_nodes.value_counts()
        .rename_axis("node")
        .reset_index(name="degree_in_export")
    )

    return node_degree.sort_values("node")


def write_cytoscape_export(
    edges: pd.DataFrame,
    filter_name: str,
    export_name: str,
) -> dict:
    """Write Cytoscape edge and node tables"""
    cyt_edges = add_cytoscape_columns(edges)
    node_table = make_node_table(cyt_edges)

    CYTOSCAPE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    edge_output = CYTOSCAPE_OUTPUT_DIR / f"{export_name}_cytoscape_edges.tsv"
    node_output = CYTOSCAPE_OUTPUT_DIR / f"{export_name}_cytoscape_nodes.tsv"

    cyt_edges.to_csv(edge_output, sep="\t", index=False)
    node_table.to_csv(node_output, sep="\t", index=False)

    return {
        "export_type": "cytoscape",
        "filter_name": filter_name,
        "export_name": export_name,
        "edge_file": str(edge_output.relative_to(PROJECT_ROOT)),
        "node_file": str(node_output.relative_to(PROJECT_ROOT)),
        "edges": len(cyt_edges),
        "nodes": len(node_table),
        "median_pident_max": cyt_edges["pident_max"].median(),
        "median_align_len_max": cyt_edges["align_len_max"].median(),
        "median_bitscore_max": cyt_edges["bitscore_max"].median(),
        "min_evalue": cyt_edges["evalue_min"].min(),
        "max_evalue": cyt_edges["evalue_min"].max(),
        "max_degree_in_export": node_table["degree_in_export"].max(),
    }


def write_consensus_table(edges: pd.DataFrame, filter_name: str) -> dict:
    """Write full consensus/all-edge table for downstream analysis"""
    CONSENSUS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    output_file = CONSENSUS_OUTPUT_DIR / f"{filter_name}.edges.tsv"
    edges.to_csv(output_file, sep="\t", index=False)

    all_nodes = pd.concat([edges["source"], edges["target"]], ignore_index=True)

    return {
        "export_type": "consensus_edge_table",
        "filter_name": filter_name,
        "export_name": filter_name,
        "edge_file": str(output_file.relative_to(PROJECT_ROOT)),
        "node_file": "",
        "edges": len(edges),
        "nodes": all_nodes.nunique(),
        "median_pident_max": edges["pident_max"].median(),
        "median_align_len_max": edges["align_len_max"].median(),
        "median_bitscore_max": edges["bitscore_max"].median(),
        "min_evalue": edges["evalue_min"].min(),
        "max_evalue": edges["evalue_min"].max(),
        "max_degree_in_export": "",
    }


def main() -> None:
    QC_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    integrated = load_integrated_edges()

    qc_rows = []

    for filter_name in NETWORK_FILTERS:
        print(f"\nCreating subset: {filter_name}")

        subset = filter_network(integrated, filter_name)

        print(f"  Full subset edges: {len(subset):,}")

        qc_rows.append(write_consensus_table(subset, filter_name))

        # Full Cytoscape exports are written only for consensus_3tools
        # Larger networks are too big for useful first visualization
        if filter_name == "consensus_3tools":
            export_name = f"evalue_1e-5_{filter_name}_full"
            print(f"  Writing full Cytoscape export: {export_name}")
            qc_rows.append(
                write_cytoscape_export(
                    edges=subset,
                    filter_name=filter_name,
                    export_name=export_name,
                )
            )

        # Visualization subsets by strongest edges
        for top_n in TOP_N_VALUES:
            top_subset = (
                subset
                .sort_values("bitscore_max", ascending=False)
                .head(top_n)
                .copy()
            )

            export_name = f"evalue_1e-5_{filter_name}_top{top_n}_bitscore"

            print(
                f"  Writing top-{top_n:,} Cytoscape export: "
                f"{len(top_subset):,} edges"
            )

            qc_rows.append(
                write_cytoscape_export(
                    edges=top_subset,
                    filter_name=filter_name,
                    export_name=export_name,
                )
            )

    qc_df = pd.DataFrame(qc_rows)

    qc_output = QC_OUTPUT_DIR / "consensus_and_cytoscape_export_qc.tsv"
    qc_df.to_csv(qc_output, sep="\t", index=False)

    print("\nConsensus and Cytoscape export QC:")
    print(qc_df.to_string(index=False))
    print(f"\nQC table written to: {qc_output}")


if __name__ == "__main__":
    main()