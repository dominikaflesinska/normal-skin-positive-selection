
#!/usr/bin/env python3
"""Summarise annotated mutation consequences by gene."""

from pathlib import Path
import sys

import pandas as pd


def load_table(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)

    header_idx = next(
        (
            i for i in range(len(raw))
            if "sampleID" in raw.iloc[i].astype(str).tolist()
        ),
        None,
    )
    if header_idx is None:
        raise ValueError("Could not locate the sampleID header.")

    df = raw.iloc[header_idx + 1:].copy()
    df.columns = raw.iloc[header_idx].tolist()
    return df.dropna(subset=["sampleID"]).copy()


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python scripts/05_gene_impact_summary.py <input.xlsx>"
        )

    source = Path(sys.argv[1])
    if not source.is_file():
        raise SystemExit(f"Input file not found: {source}")

    df = load_table(source)
    annotated = df.dropna(subset=["gene_name", "impact"]).copy()

    summary = pd.crosstab(
        annotated["gene_name"],
        annotated["impact"],
    ).sort_index()

    summary["Total annotated records"] = summary.sum(axis=1)
    summary = summary.sort_values(
        "Total annotated records", ascending=False
    )

    results = Path("results")
    results.mkdir(parents=True, exist_ok=True)

    summary.to_csv(results / "gene_impact_counts.csv")

    report = f"""# Gene-level mutation impact summary

## Dataset

- Total source records: {len(df)}
- Records with both gene and impact annotations: {len(annotated)}
- Genes represented among annotated records: {annotated['gene_name'].nunique()}
- Records excluded because gene or impact annotation was missing: {len(df) - len(annotated)}

## Genes with the most annotated records

{summary.head(15).to_string()}

## Interpretation limits

These are counts of previously called variants in the study's targeted dataset.
They are not normalized for gene length, sequence context, mutation opportunity,
or sample-specific effects. The counts do not constitute a dN/dS estimate or
a statistical test of positive selection.

Genes with more records should not automatically be interpreted as more strongly
selected.
"""

    (results / "gene_impact_summary.md").write_text(
        report, encoding="utf-8"
    )

    print(report)
    print(f"Wrote {results / 'gene_impact_counts.csv'}")
    print(f"Wrote {results / 'gene_impact_summary.md'}")


if __name__ == "__main__":
    main()
