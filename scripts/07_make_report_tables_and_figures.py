"""
Create report-ready tables and figures from QC outputs and consensus exports.

Usage:
    python scripts/07_make_report_tables_and_figures.py

Required input files:
    - data/quality_checks/edge_table_construction_qc.tsv
    - data/quality_checks/integrated_tool_support_qc.tsv
    - data/quality_checks/integrated_network_summary_statistics.tsv
    - data/quality_checks/consensus_and_cytoscape_export_qc.tsv
    - data/edges/consensus/best_bitscore/evalue_1e-5/*.edges.tsv

All of these files are produced by the earlier scripts in the workflow. If any of
them are missing, this final reporting step will not be able to build the summary
tables and figures.

This script reads QC tables and consensus edge exports, produces the final report
tables and individual figures, and combines selected plots into the numbered
multi-panel figures used in the assignment report.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import math
import shutil


PROJECT_ROOT = Path(__file__).resolve().parents[1]

QC_DIR = PROJECT_ROOT / "data" / "quality_checks"
RESULTS_DIR = PROJECT_ROOT / "results"
TABLE_DIR = RESULTS_DIR / "tables"
FIGURE_DIR = RESULTS_DIR / "figures"

CONSENSUS_EDGE_DIR = (
    PROJECT_ROOT
    / "data"
    / "edges"
    / "consensus"
    / "best_bitscore"
    / "evalue_1e-5"
)

THRESHOLD_ORDER = ["raw", "evalue_1e-3", "evalue_1e-5", "evalue_1e-10"]
NETWORK_FILTER_ORDER = ["all_edges", "consensus_2plus", "consensus_3tools"]

MAGENTA = "#ff00cc"
CYAN = "#00d5ff"
LIGHT_GREY = "#cfcfcf"
DARK_GREY = "#555555"

TOOL_COLORS = {
    "blastp": MAGENTA,
    "mmseqs": CYAN,
    "swipe": LIGHT_GREY,
}

SUPPORT_COLORS = {
    "edges_supported_by_1_tool": CYAN,
    "edges_supported_by_2_tools": LIGHT_GREY,
    "edges_supported_by_3_tools": MAGENTA,
}

NETWORK_FILTER_COLORS = {
    "all_edges": CYAN,
    "consensus_2plus": LIGHT_GREY,
    "consensus_3tools": MAGENTA,
}

QUALITY_COLORS = {
    "median_pident_max": MAGENTA,
    "median_align_len_max": CYAN,
    "median_bitscore_max": LIGHT_GREY,
}


class UnionFind:
    """Memory-friendly connected component calculation for undirected graphs"""

    def __init__(self):
        self.parent = {}
        self.size = {}

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 1

    def find(self, x: str) -> str:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        self.add(a)
        self.add(b)

        root_a = self.find(a)
        root_b = self.find(b)

        if root_a == root_b:
            return

        if self.size[root_a] < self.size[root_b]:
            root_a, root_b = root_b, root_a

        self.parent[root_b] = root_a
        self.size[root_a] += self.size[root_b]

    def component_sizes(self) -> list[int]:
        counts = {}

        for node in self.parent:
            root = self.find(node)
            counts[root] = counts.get(root, 0) + 1

        return list(counts.values())


def setup_output_dirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def save_table(df: pd.DataFrame, filename: str) -> None:
    output_path = TABLE_DIR / filename
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Table written: {output_path}")


def save_figure(filename: str) -> None:
    output_path = FIGURE_DIR / filename
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Figure written: {output_path}")
    plt.close()


def prepare_edge_count_table() -> pd.DataFrame:
    qc = pd.read_csv(QC_DIR / "edge_table_construction_qc.tsv", sep="\t")

    df = qc[
        qc["collapse_rule"] == "best_bitscore"
    ][
        [
            "tool",
            "threshold",
            "alignment_rows_after_threshold",
            "unique_edges",
            "collapsed_duplicate_hits",
            "median_pident",
            "median_length",
            "median_bitscore",
        ]
    ].copy()

    df["threshold"] = pd.Categorical(
        df["threshold"],
        categories=THRESHOLD_ORDER,
        ordered=True,
    )

    df = df.sort_values(["threshold", "tool"])

    save_table(df, "edge_counts_by_tool_threshold.tsv")
    return df


def plot_edge_counts_by_tool(edge_counts: pd.DataFrame) -> None:
    pivot = edge_counts.pivot(
        index="threshold",
        columns="tool",
        values="unique_edges",
    ).loc[THRESHOLD_ORDER]

    colors = [TOOL_COLORS[col] for col in pivot.columns]

    ax = pivot.plot(kind="bar", figsize=(9, 5), color=colors)

    ax.set_title("Unique protein-pair edges by tool and E-value threshold")
    ax.set_xlabel("E-value threshold")
    ax.set_ylabel("Unique edges")
    ax.tick_params(axis="x", rotation=30)

    save_figure("edge_counts_by_tool_threshold.png")


def prepare_integrated_tool_support_table() -> pd.DataFrame:
    qc = pd.read_csv(QC_DIR / "integrated_tool_support_qc.tsv", sep="\t")

    df = qc[
        qc["collapse_rule"] == "best_bitscore"
    ][
        [
            "threshold",
            "integrated_edges",
            "edges_supported_by_1_tool",
            "edges_supported_by_2_tools",
            "edges_supported_by_3_tools",
            "blastp_only",
            "mmseqs_only",
            "swipe_only",
            "blastp_mmseqs",
            "blastp_swipe",
            "mmseqs_swipe",
            "blastp_mmseqs_swipe",
            "median_pident_max",
            "median_align_len_max",
            "median_bitscore_max",
        ]
    ].copy()

    df["threshold"] = pd.Categorical(
        df["threshold"],
        categories=THRESHOLD_ORDER,
        ordered=True,
    )

    df = df.sort_values("threshold")

    df["fraction_supported_by_1_tool"] = (
        df["edges_supported_by_1_tool"] / df["integrated_edges"]
    )
    df["fraction_supported_by_2_tools"] = (
        df["edges_supported_by_2_tools"] / df["integrated_edges"]
    )
    df["fraction_supported_by_3_tools"] = (
        df["edges_supported_by_3_tools"] / df["integrated_edges"]
    )

    save_table(df, "integrated_tool_support_summary.tsv")
    return df


def plot_integrated_tool_support_absolute(tool_support: pd.DataFrame) -> None:
    plot_df = tool_support.set_index("threshold")[
        [
            "edges_supported_by_1_tool",
            "edges_supported_by_2_tools",
            "edges_supported_by_3_tools",
        ]
    ].loc[THRESHOLD_ORDER]

    colors = [SUPPORT_COLORS[col] for col in plot_df.columns]

    ax = plot_df.plot(kind="bar", stacked=True, figsize=(9, 5), color=colors)

    ax.set_title("Cross-tool support of integrated protein similarity edges")
    ax.set_xlabel("E-value threshold")
    ax.set_ylabel("Integrated edges")
    ax.tick_params(axis="x", rotation=30)

    save_figure("integrated_tool_support_by_threshold.png")


def plot_integrated_tool_support_fraction(tool_support: pd.DataFrame) -> None:
    plot_df = tool_support.set_index("threshold")[
        [
            "fraction_supported_by_1_tool",
            "fraction_supported_by_2_tools",
            "fraction_supported_by_3_tools",
        ]
    ].loc[THRESHOLD_ORDER]

    renamed = plot_df.rename(
        columns={
            "fraction_supported_by_1_tool": "1 tool",
            "fraction_supported_by_2_tools": "2 tools",
            "fraction_supported_by_3_tools": "3 tools",
        }
    )

    colors = [CYAN, LIGHT_GREY, MAGENTA]

    ax = renamed.plot(kind="bar", stacked=True, figsize=(9, 5), color=colors)

    ax.set_title("Fractional cross-tool support by E-value threshold")
    ax.set_xlabel("E-value threshold")
    ax.set_ylabel("Fraction of integrated edges")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="x", rotation=30)

    save_figure("integrated_tool_support_fraction_by_threshold.png")


def plot_support_pattern_primary_threshold(tool_support: pd.DataFrame) -> None:
    row = tool_support[tool_support["threshold"] == "evalue_1e-5"].iloc[0]

    pattern_counts = pd.Series(
        {
            "BLASTP only": row["blastp_only"],
            "MMseqs2 only": row["mmseqs_only"],
            "SWIPE only": row["swipe_only"],
            "BLASTP + MMseqs2": row["blastp_mmseqs"],
            "BLASTP + SWIPE": row["blastp_swipe"],
            "MMseqs2 + SWIPE": row["mmseqs_swipe"],
            "All three tools": row["blastp_mmseqs_swipe"],
        }
    ).sort_values(ascending=True)

    colors = []
    for label in pattern_counts.index:
        if label == "All three tools":
            colors.append(MAGENTA)
        elif "only" in label:
            colors.append(CYAN)
        else:
            colors.append(LIGHT_GREY)

    ax = pattern_counts.plot(kind="barh", figsize=(9, 5), color=colors)

    ax.set_title("Support-pattern composition at E-value <= 1e-5")
    ax.set_xlabel("Integrated edges")
    ax.set_ylabel("Support pattern")

    save_figure("support_pattern_composition_evalue_1e-5.png")


def prepare_primary_integrated_network_summary() -> pd.DataFrame:
    qc = pd.read_csv(QC_DIR / "integrated_network_summary_statistics.tsv", sep="\t")

    df = qc[
        (qc["collapse_rule"] == "best_bitscore")
        & (qc["threshold"] == "evalue_1e-5")
    ][
        [
            "threshold",
            "network_filter",
            "nodes_with_at_least_one_edge",
            "edges",
            "average_degree",
            "connected_components",
            "largest_component_size",
            "largest_component_fraction",
            "median_degree",
            "max_degree",
            "median_pident_max",
            "median_align_len_max",
            "median_bitscore_max",
        ]
    ].copy()

    df["network_filter"] = pd.Categorical(
        df["network_filter"],
        categories=NETWORK_FILTER_ORDER,
        ordered=True,
    )

    df = df.sort_values("network_filter")

    save_table(df, "primary_integrated_network_summary.tsv")
    return df


def plot_largest_component_by_filter(primary_summary: pd.DataFrame) -> None:
    plot_df = primary_summary.set_index("network_filter").loc[NETWORK_FILTER_ORDER]

    colors = [NETWORK_FILTER_COLORS[idx] for idx in plot_df.index]

    ax = plot_df["largest_component_fraction"].plot(
        kind="bar",
        figsize=(8, 5),
        color=colors,
    )

    ax.set_title("Largest component fraction by consensus filter")
    ax.set_xlabel("Network filter")
    ax.set_ylabel("Largest component fraction")
    ax.tick_params(axis="x", rotation=30)

    save_figure("largest_component_by_network_filter.png")


def plot_alignment_quality_by_filter(primary_summary: pd.DataFrame) -> None:
    plot_df = primary_summary.set_index("network_filter").loc[NETWORK_FILTER_ORDER]

    quality_df = plot_df[
        [
            "median_pident_max",
            "median_align_len_max",
            "median_bitscore_max",
        ]
    ]

    colors = [MAGENTA, CYAN, LIGHT_GREY]

    ax = quality_df.plot(kind="bar", figsize=(9, 5), color=colors)

    ax.set_title("Alignment quality metrics by consensus filter")
    ax.set_xlabel("Network filter")
    ax.set_ylabel("Median value")
    ax.tick_params(axis="x", rotation=30)

    save_figure("alignment_quality_by_consensus_filter.png")


def plot_relative_alignment_quality_by_filter(primary_summary: pd.DataFrame) -> None:
    plot_df = primary_summary.set_index("network_filter").loc[NETWORK_FILTER_ORDER]

    quality_df = plot_df[
        [
            "median_pident_max",
            "median_align_len_max",
            "median_bitscore_max",
        ]
    ].copy()

    relative_df = quality_df / quality_df.loc["all_edges"]

    renamed = relative_df.rename(
        columns={
            "median_pident_max": "median pident",
            "median_align_len_max": "median alignment length",
            "median_bitscore_max": "median bitscore",
        }
    )

    colors = [MAGENTA, CYAN, LIGHT_GREY]

    ax = renamed.plot(kind="bar", figsize=(9, 5), color=colors)

    ax.axhline(1.0, color=DARK_GREY, linewidth=1, linestyle="--")

    ax.set_title("Relative alignment-quality increase by consensus filter")
    ax.set_xlabel("Network filter")
    ax.set_ylabel("Relative value compared with all_edges")
    ax.tick_params(axis="x", rotation=30)

    save_figure("relative_alignment_quality_by_consensus_filter.png")


def load_consensus_edges(filter_name: str) -> pd.DataFrame:
    input_file = CONSENSUS_EDGE_DIR / f"{filter_name}.edges.tsv"

    if not input_file.exists():
        raise FileNotFoundError(f"Missing consensus edge file: {input_file}")

    usecols = [
        "source",
        "target",
        "n_tools",
        "support_pattern",
        "pident_max",
        "align_len_max",
        "evalue_min",
        "bitscore_max",
    ]

    return pd.read_csv(input_file, sep="\t", usecols=usecols)


def calculate_degree_series(edges: pd.DataFrame) -> pd.Series:
    all_nodes = pd.concat([edges["source"], edges["target"]], ignore_index=True)
    return all_nodes.value_counts()


def make_degree_distribution_table() -> pd.DataFrame:
    rows = []

    for filter_name in NETWORK_FILTER_ORDER:
        print(f"Calculating degree distribution: {filter_name}")
        edges = load_consensus_edges(filter_name)
        degree_series = calculate_degree_series(edges)

        degree_counts = degree_series.value_counts().sort_index()

        for degree, node_count in degree_counts.items():
            rows.append(
                {
                    "network_filter": filter_name,
                    "degree": int(degree),
                    "node_count": int(node_count),
                }
            )

    df = pd.DataFrame(rows)
    save_table(df, "degree_distribution_primary_networks.tsv")
    return df


def plot_degree_ccdf(degree_distribution: pd.DataFrame) -> None:
    plt.figure(figsize=(8, 5))

    for filter_name in NETWORK_FILTER_ORDER:
        sub = degree_distribution[
            degree_distribution["network_filter"] == filter_name
        ].sort_values("degree")

        total_nodes = sub["node_count"].sum()

        # CCDF: fraction of nodes with degree >= k
        sub = sub.copy()
        sub["nodes_with_degree_at_least_k"] = sub["node_count"][::-1].cumsum()[::-1]
        sub["fraction_nodes_with_degree_at_least_k"] = (
            sub["nodes_with_degree_at_least_k"] / total_nodes
        )

        plt.plot(
            sub["degree"],
            sub["fraction_nodes_with_degree_at_least_k"],
            marker="o",
            markersize=3,
            linewidth=1.5,
            label=filter_name,
            color=NETWORK_FILTER_COLORS[filter_name],
        )

    plt.xscale("log")
    plt.yscale("log")
    plt.title("Degree distribution CCDF of primary integrated networks")
    plt.xlabel("Degree k")
    plt.ylabel("Fraction of nodes with degree ≥ k")
    plt.legend()

    save_figure("degree_ccdf_primary_integrated_networks.png")


def calculate_component_sizes(edges: pd.DataFrame) -> list[int]:
    uf = UnionFind()

    for source, target in zip(edges["source"], edges["target"]):
        uf.union(source, target)

    return uf.component_sizes()


def make_component_size_distribution_table() -> pd.DataFrame:
    rows = []

    for filter_name in NETWORK_FILTER_ORDER:
        print(f"Calculating component size distribution: {filter_name}")
        edges = load_consensus_edges(filter_name)
        component_sizes = calculate_component_sizes(edges)

        size_counts = pd.Series(component_sizes).value_counts().sort_index()

        for component_size, component_count in size_counts.items():
            rows.append(
                {
                    "network_filter": filter_name,
                    "component_size": int(component_size),
                    "component_count": int(component_count),
                }
            )

    df = pd.DataFrame(rows)
    save_table(df, "component_size_distribution_primary_networks.tsv")
    return df


def plot_component_size_distribution(component_distribution: pd.DataFrame) -> None:
    plt.figure(figsize=(8, 5))

    for filter_name in NETWORK_FILTER_ORDER:
        sub = component_distribution[
            component_distribution["network_filter"] == filter_name
        ].sort_values("component_size")

        plt.scatter(
            sub["component_size"],
            sub["component_count"],
            s=18,
            label=filter_name,
            color=NETWORK_FILTER_COLORS[filter_name],
            alpha=0.8,
        )

    plt.xscale("log")
    plt.yscale("log")
    plt.title("Component size distribution of primary integrated networks")
    plt.xlabel("Component size")
    plt.ylabel("Number of components")
    plt.legend()

    save_figure("component_size_distribution_primary_networks.png")

def make_panel_figure(
    image_filenames: list[str],
    panel_titles: list[str],
    output_filename: str,
    ncols: int = 2,
    figsize: tuple[int, int] = (14, 8),
) -> None:
    """Combine existing PNG figures into one multi-panel report figure"""
    n_panels = len(image_filenames)
    nrows = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

    if nrows == 1 and ncols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, (filename, title) in enumerate(zip(image_filenames, panel_titles)):
        image_path = FIGURE_DIR / filename

        if not image_path.exists():
            raise FileNotFoundError(f"Missing figure for panel: {image_path}")

        img = plt.imread(image_path)
        axes[i].imshow(img)
        axes[i].axis("off")
        axes[i].set_title(
            f"{chr(65 + i)}. {title}",
            loc="left",
            fontsize=12,
            fontweight="bold",
        )

    for j in range(n_panels, len(axes)):
        axes[j].axis("off")

    output_path = FIGURE_DIR / output_filename
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Panel figure written: {output_path}")

def copy_report_figure(source_filename: str, output_filename: str) -> None:
    """Copy an existing single-panel figure to a report-numbered filename."""
    source_path = FIGURE_DIR / source_filename
    output_path = FIGURE_DIR / output_filename

    if not source_path.exists():
        raise FileNotFoundError(f"Missing source figure: {source_path}")

    shutil.copyfile(source_path, output_path)
    print(f"Single-panel report figure written: {output_path}")

def make_combined_report_figures() -> None:
    """Create report-ready figures from individual plots."""

    make_panel_figure(
        image_filenames=[
            "edge_counts_by_tool_threshold.png",
            "integrated_tool_support_fraction_by_threshold.png",
        ],
        panel_titles=[
            "Unique edges by tool and threshold",
            "Fractional cross-tool support",
        ],
        output_filename="figure_1_tool_threshold_and_support.png",
        ncols=2,
        figsize=(16, 7),
    )

    make_panel_figure(
        image_filenames=[
            "support_pattern_composition_evalue_1e-5.png",
            "relative_alignment_quality_by_consensus_filter.png",
        ],
        panel_titles=[
            "Support-pattern composition at E-value ≤ 1e-5",
            "Relative alignment-quality increase",
        ],
        output_filename="figure_2_primary_support_and_quality.png",
        ncols=2,
        figsize=(16, 7),
    )

    copy_report_figure(
        source_filename="degree_ccdf_primary_integrated_networks.png",
        output_filename="figure_3_degree_distribution_ccdf.png",
    )

    make_panel_figure(
        image_filenames=[
            "largest_component_by_network_filter.png",
            "component_size_distribution_primary_networks.png",
        ],
        panel_titles=[
            "Largest component fraction",
            "Component size distribution",
        ],
        output_filename="figure_4_network_fragmentation_and_components.png",
        ncols=1,
        figsize=(9, 12),
    )

def main() -> None:
    setup_output_dirs()

    edge_counts = prepare_edge_count_table()
    plot_edge_counts_by_tool(edge_counts)

    tool_support = prepare_integrated_tool_support_table()
    plot_integrated_tool_support_absolute(tool_support)
    plot_integrated_tool_support_fraction(tool_support)
    plot_support_pattern_primary_threshold(tool_support)

    primary_summary = prepare_primary_integrated_network_summary()
    plot_largest_component_by_filter(primary_summary)
    plot_alignment_quality_by_filter(primary_summary)
    plot_relative_alignment_quality_by_filter(primary_summary)

    degree_distribution = make_degree_distribution_table()
    plot_degree_ccdf(degree_distribution)

    component_distribution = make_component_size_distribution_table()
    plot_component_size_distribution(component_distribution)

    make_combined_report_figures()

    print("\nReport tables and figures completed.")

if __name__ == "__main__":
    main()