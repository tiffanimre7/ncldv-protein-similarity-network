# Comparative construction and analysis of NCLDV protein similarity networks

This repository contains the analysis workflow for constructing and comparing protein similarity networks from all-vs-all BLASTP, MMseqs2 and SWIPE searches of Nucleocytoviricota/NCLDV protein sequences.

The project evaluates how sequence comparison tool choice, E-value thresholding, bit-score-based edge selection and cross-tool consensus influence the resulting protein similarity networks.

## Project overview

Large DNA viruses have complex and variable gene repertoires. Because not all viral proteins share universal marker genes, network-based approaches can complement traditional phylogenetic analyses by representing many local sequence similarity relationships.

In this project, protein similarity networks were constructed from standardized all-vs-all alignment outputs. Nodes represent proteins and edges represent detected sequence similarity relationships between protein pairs. Edges were compared across BLASTP, MMseqs2 and SWIPE, filtered under multiple E-value thresholds and integrated into support-based consensus networks.

## Main workflow

1. Standardize BLASTP, MMseqs2 and SWIPE alignment outputs.
2. Remove self-hits and normalize sequence identifiers.
3. Convert alignments into undirected protein-pair edge tables.
4. Collapse duplicate protein-pair hits primarily by highest bit score.
5. Generate networks under multiple E-value thresholds.
6. Integrate tool support across BLASTP, MMseqs2 and SWIPE.
7. Analyze edge counts, cross-tool support, degree distributions and connected components.
8. Export report-ready figures and Cytoscape visualization files.

## Main result

The three tools differed strongly in the number of recovered similarity edges, especially in the raw outputs. Stricter E-value filtering reduced these differences and increased cross-tool agreement. At the primary threshold of E-value <= 1e-5, the integrated all-edges network contained 1,223,039 edges, the consensus_2plus network contained 999,485 edges, and the consensus_3tools network contained 735,449 edges.

## What is included

- `scripts/`: reproducible analysis scripts.
- `results/tables/`: report-ready summary tables.
- `results/figures/`: final figures and supporting plots.
- `docs/report.pdf`: the assignment report.
- `data/README.md`: data provenance and exclusion note.
- `environment/requirements.txt`: minimal Python dependencies.

## What is intentionally excluded

This repository does not include the raw all-vs-all alignment outputs or the large intermediate edge tables generated from the course/HPC environment.

## Reproducibility

The workflow is implemented as a sequence of numbered scripts:

- `scripts/01_standardize_alignments.py`
- `scripts/02_build_edge_tables.py`
- `scripts/03_network_summary_statistics.py`
- `scripts/04_integrate_tool_support.py`
- `scripts/05_integrated_network_summary.py`
- `scripts/06_export_consensus_and_cytoscape.py`
- `scripts/07_make_report_tables_and_figures.py`

The core Python dependencies are listed in `environment/requirements.txt`.

## Reproduce Workflow

```bash
python -m venv .venv
pip install -r environment/requirements.txt
python scripts/01_standardize_alignments.py
python scripts/02_build_edge_tables.py
python scripts/03_network_summary_statistics.py
python scripts/04_integrate_tool_support.py
python scripts/05_integrated_network_summary.py
python scripts/06_export_consensus_and_cytoscape.py
python scripts/07_make_report_tables_and_figures.py
```

### Upstream HPC alignment step

The upstream all-vs-all alignment step was run on the LiSC cluster using `scripts/00_run_hpc_alignments.sh`.

This script executes BLASTP, SWIPE and MMseqs2, removes self-hits, records runtime statistics and writes the final `.m8` files used for downstream analysis. Large raw alignment outputs are not included in this repository because of file size constraints.

## Suggested repository layout

```text
ncldv-protein-similarity-network/
├── README.md
├── scripts/
├── results/
│   ├── figures/
│   └── tables/
├── docs/
│   └── report.pdf
├── data/
│   └── README.md
├── environment/
│   └── requirements.txt
└── .gitignore
```
