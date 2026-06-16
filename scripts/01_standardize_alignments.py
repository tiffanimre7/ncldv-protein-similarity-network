"""
Standardize BLASTP, MMseqs2, and SWIPE alignment outputs and write a QC summary.

Usage:
    python scripts/01_standardize_alignments.py

Required input files:
    - data/ncldv.blastp.wo_self.m8
    - data/ncldv.mmseqs.wo_self.m8
    - data/ncldv.swipe.cleaned.wo_self.m8

If you are using your own data, place the raw alignment outputs in the repository's
data/ directory with the same file names, or update the INPUTS paths in this script.

This script standardizes the output of the three alignment tools to a common format,
removes self-hits, performs basic QC checks, and saves the standardized alignments
plus a summary TSV under data/quality_checks/.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUTS = {
    "blastp": PROJECT_ROOT / "data" / "ncldv.blastp.wo_self.m8",
    "mmseqs": PROJECT_ROOT / "data" / "ncldv.mmseqs.wo_self.m8",
    "swipe": PROJECT_ROOT / "data" / "ncldv.swipe.cleaned.wo_self.m8",
}

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed"
QC_OUTPUT = PROJECT_ROOT / "data" / "quality_checks"

def clean_id(value: str) -> str:
    """Clean sequence IDs by removing prefixes and suffixes"""
    value = str(value)
    value = value.replace("lcl|", "")
    if value.startswith("ref|") and value.endswith("|"):
        value = value[4:-1]
    return value


def standardize_alignment(tool: str, input_path: Path, output_path: Path) -> dict:
    """Standardize alignment output to a common format and perform QC checks"""
    print(f"Processing {tool}: {input_path}") # Debug statement to indicate which tool and file is being processed
    
    df = pd.read_csv(
        input_path,
        sep="\t",
        names=COLUMNS,
        header=None,
        dtype={
            "qseqid": str,
            "sseqid": str,
        },
    )

    raw_rows = len(df)

    df["qseqid"] = df["qseqid"].map(clean_id)
    df["sseqid"] = df["sseqid"].map(clean_id)

    numeric_cols = [
        "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    missing_numeric_rows = df[numeric_cols].isna().any(axis=1).sum()

    if df["pident"].max() <= 1.0:
        df["pident"] = df["pident"] * 100

    self_hits_before_filter = (df["qseqid"] == df["sseqid"]).sum()
    df = df[df["qseqid"] != df["sseqid"]].copy()

    df.insert(0, "tool", tool)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    print(f"Finished {tool}: {len(df):,} standardized rows written to {output_path}") # Debug statement to indicate completion of processing for the tool and number of rows written

    return {
        "tool": tool,
        "input_file": str(input_path),
        "output_file": str(output_path),
        "raw_rows": raw_rows,
        "self_hits_removed": int(self_hits_before_filter),
        "rows_after_filtering": len(df),
        "missing_numeric_rows": int(missing_numeric_rows),
        "min_pident": df["pident"].min(),
        "max_pident": df["pident"].max(),
        "min_evalue": df["evalue"].min(),
        "max_evalue": df["evalue"].max(),
        "min_bitscore": df["bitscore"].min(),
        "max_bitscore": df["bitscore"].max(),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_OUTPUT.mkdir(parents=True, exist_ok=True)

    qc_rows = []

    for tool, input_path in INPUTS.items():
        output_path = OUTPUT_DIR / f"{tool}.standardized.tsv"
        qc_rows.append(standardize_alignment(tool, input_path, output_path))

    qc_df = pd.DataFrame(qc_rows)

    qc_output = QC_OUTPUT / "alignment_standardization_qc.tsv"
    qc_df.to_csv(qc_output, sep="\t", index=False)

    print("\nAlignment standardization QC:")
    print(qc_df.to_string(index=False))
    print(f"\nQC table written to: {qc_output}")


if __name__ == "__main__":
    main()