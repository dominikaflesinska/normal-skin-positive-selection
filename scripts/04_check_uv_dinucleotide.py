
#!/usr/bin/env python3
"""Inspect multi-nucleotide changes relevant to UV-associated mutagenesis."""

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
            "Usage: python scripts/04_check_uv_dinucleotide.py <input.xlsx>"
        )

    source = Path(sys.argv[1])
    if not source.is_file():
        raise SystemExit(f"Input file not found: {source}")

    df = load_table(source)

    ref = df["ref_nt"].astype(str).str.upper().str.strip()
    alt = df["mut_nt"].astype(str).str.upper().str.strip()

    multi = (
        ref.str.len().gt(1)
        & alt.str.len().eq(ref.str.len())
        & ref.str.fullmatch("[ACGT]+", na=False)
        & alt.str.fullmatch("[ACGT]+", na=False)
        & ref.ne(alt)
    )

    changes = df.loc[multi].copy()
    changes["allele_change"] = (
        ref.loc[multi] + ">" + alt.loc[multi]
    )

    change_counts = (
        changes["allele_change"]
        .value_counts()
        .rename_axis("allele_change")
        .rename("count")
    )

    cc_tt = changes.loc[
        changes["allele_change"].eq("CC>TT")
    ]

    results = Path("results")
    results.mkdir(parents=True, exist_ok=True)

    change_counts.to_csv(results / "multi_nucleotide_change_counts.csv")

    report = f"""# Multi-nucleotide change inspection

## Dataset

- Total source records: {len(df)}
- Equal-length multi-nucleotide changes with A/C/G/T alleles: {int(multi.sum())}
- Records represented as `CC>TT`: {len(cc_tt)}

## Most frequent allele changes

{change_counts.head(20).to_string()}

## Interpretation limits

Allele strings alone do not establish that a change is a validated dinucleotide
mutation or that it was caused by ultraviolet radiation. The source table's
variant representation, genomic coordinates, and reference sequence should be
checked before interpreting adjacent substitutions as a single event.

This report is descriptive and does not test a UV mutational signature.
"""

    (results / "uv_dinucleotide_check.md").write_text(
        report, encoding="utf-8"
    )

    print(report)
    print(f"Wrote {results / 'multi_nucleotide_change_counts.csv'}")
    print(f"Wrote {results / 'uv_dinucleotide_check.md'}")


if __name__ == "__main__":
    main()
