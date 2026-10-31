
#!/usr/bin/env python3
"""Classify variant records by reference and alternate allele length."""

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


def classify_variant(ref: object, alt: object) -> str:
    if pd.isna(ref) or pd.isna(alt):
        return "Missing allele"

    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()

    if not ref or not alt:
        return "Missing allele"

    if ref == "-" and alt == "-":
        return "Both alleles are gaps"

    if ref == "-":
        return "Insertion-like representation"

    if alt == "-":
        return "Deletion-like representation"

    if ref == alt:
        return "Identical reference and alternate alleles"

    if len(ref) == 1 and len(alt) == 1:
        if ref in "ACGT" and alt in "ACGT":
            return "Single-nucleotide substitution"

    if len(ref) == len(alt) and set(ref + alt) <= set("ACGT"):
        return "Equal-length multi-nucleotide change"

    if set(ref + alt) <= set("ACGT"):
        return "Unequal-length sequence change"

    return "Other or ambiguous representation"


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python scripts/03_classify_variants.py <input.xlsx>"
        )

    source = Path(sys.argv[1])
    if not source.is_file():
        raise SystemExit(f"Input file not found: {source}")

    df = load_table(source)

    df["variant_class"] = [
        classify_variant(ref, alt)
        for ref, alt in zip(df["ref_nt"], df["mut_nt"])
    ]

    counts = (
        df["variant_class"]
        .value_counts()
        .rename_axis("variant_class")
        .rename("count")
    )

    results = Path("results")
    results.mkdir(parents=True, exist_ok=True)

    counts.to_csv(results / "variant_class_counts.csv")

    report = f"""# Variant classification

## Dataset

- Source file: `{source.name}`
- Total records: {len(df)}
- Number of classified records: {int(counts.sum())}
- Classification totals match the input record count: {int(counts.sum()) == len(df)}

## Variant classes

{counts.to_string()}

## Method and limitations

Records were classified using the lengths and literal contents of the
reference and alternate allele fields. Gap symbols (`-`) were treated as
insertion-like or deletion-like representations.

These labels describe the table representation only. They do not replace
normalization against a reference genome or validation of individual
variants. Equal-length multi-nucleotide changes are not necessarily
independent single-nucleotide events.

The source table contains previously called variants. This script does not
re-call variants, infer mutational mechanisms, or estimate positive selection.
"""

    (results / "variant_classification.md").write_text(
        report, encoding="utf-8"
    )

    print(report)
    print(f"Wrote {results / 'variant_class_counts.csv'}")
    print(f"Wrote {results / 'variant_classification.md'}")


if __name__ == "__main__":
    main()
