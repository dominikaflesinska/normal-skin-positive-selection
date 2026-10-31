
#!/usr/bin/env python3
"""Summarise single-nucleotide substitutions in normal human skin."""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BASES = set("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SUBSTITUTIONS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


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
            "Usage: python scripts/02_mutational_spectrum.py <input.xlsx>"
        )

    source = Path(sys.argv[1])
    if not source.is_file():
        raise SystemExit(f"Input file not found: {source}")

    df = load_table(source)

    ref = df["ref_nt"].astype(str).str.upper()
    alt = df["mut_nt"].astype(str).str.upper()
    context = df["trinuc"].astype(str).str.upper()

    is_snv = (
        ref.str.len().eq(1)
        & alt.str.len().eq(1)
        & ref.isin(BASES)
        & alt.isin(BASES)
        & ref.ne(alt)
    )

    snv = df.loc[is_snv].copy()
    ref = ref.loc[is_snv]
    alt = alt.loc[is_snv]
    context = context.loc[is_snv]

    valid_context = (
        context.str.fullmatch("[ACGT]{3}", na=False)
        & context.str[1].eq(ref)
    )

    snv = snv.loc[valid_context].copy()
    ref = ref.loc[valid_context]
    alt = alt.loc[valid_context]
    context = context.loc[valid_context]

    normalized_context = []
    normalized_substitution = []

    for reference, alternate, trinucleotide in zip(ref, alt, context):
        if reference in {"A", "G"}:
            trinucleotide = reverse_complement(trinucleotide)
            reference = reference.translate(COMPLEMENT)
            alternate = alternate.translate(COMPLEMENT)

        normalized_context.append(trinucleotide)
        normalized_substitution.append(f"{reference}>{alternate}")

    snv["normalized_context"] = normalized_context
    snv["substitution"] = normalized_substitution
    snv["channel"] = (
        snv["normalized_context"].str[0]
        + "["
        + snv["substitution"]
        + "]"
        + snv["normalized_context"].str[2]
    )

    contexts = [
        left + "[" + substitution + "]" + right
        for substitution in SUBSTITUTIONS
        for left in "ACGT"
        for right in "ACGT"
    ]

    channel_counts = (
        snv["channel"].value_counts()
        .reindex(contexts, fill_value=0)
        .rename_axis("mutation_context")
        .rename("count")
    )

    substitution_counts = (
        snv["substitution"].value_counts()
        .reindex(SUBSTITUTIONS, fill_value=0)
        .rename_axis("substitution")
        .rename("count")
    )

    results = Path("results")
    results.mkdir(parents=True, exist_ok=True)

    channel_counts.to_csv(results / "96_context_counts.csv")
    substitution_counts.to_csv(results / "six_class_counts.csv")

    fig, ax = plt.subplots(figsize=(16, 5))
    channel_counts.plot(kind="bar", ax=ax, width=0.9)
    ax.set_xlabel("Normalized trinucleotide context")
    ax.set_ylabel("Number of SNVs")
    ax.set_title("Single-nucleotide substitution spectrum (96 contexts)")
    ax.tick_params(axis="x", labelsize=6, rotation=90)
    fig.tight_layout()
    fig.savefig(results / "mutational_spectrum_96.png", dpi=200)
    plt.close(fig)

    report = f"""# Mutational spectrum

## Dataset and quality checks

- Source file: `{source.name}`
- Records in source table: {len(df)}
- Single-nucleotide substitutions: {int(is_snv.sum())}
- SNVs with valid three-base context matching the reference allele: {len(snv)}
- SNVs excluded from context analysis: {int(is_snv.sum()) - len(snv)}

## Six substitution classes

{substitution_counts.to_string()}

## Trinucleotide contexts

The 96 channels use pyrimidine-oriented substitutions (C or T as the
reference base), with the adjacent reference bases retained on both sides.
Counts are descriptive and have not been normalized for sequence opportunity
or compared with a formal mutational-signature reference.

## Interpretation limits

The source is a filtered set of previously called mutations in targeted genes,
not a genome-wide catalogue of mutations. These counts alone do not establish
a UV mutational signature and do not estimate positive selection.
"""

    (results / "mutational_spectrum.md").write_text(
        report, encoding="utf-8"
    )

    print(report)
    print(f"Wrote {results / 'six_class_counts.csv'}")
    print(f"Wrote {results / '96_context_counts.csv'}")
    print(f"Wrote {results / 'mutational_spectrum_96.png'}")
    print(f"Wrote {results / 'mutational_spectrum.md'}")


if __name__ == "__main__":
    main()

