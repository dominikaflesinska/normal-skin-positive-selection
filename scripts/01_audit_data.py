
#!/usr/bin/env python3
"""Structural audit of the Martincorena et al. (2015) mutation table."""

from pathlib import Path
import sys

import pandas as pd


REQUIRED_COLUMNS = {
    "sampleID", "chr", "pos", "ref_nt", "mut_nt",
    "reads", "cov", "vaf", "trinuc", "qval",
    "gene_name", "impact",
}


def load_mutations(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)

    header_idx = next(
        (
            i for i in range(len(raw))
            if "sampleID" in raw.iloc[i].astype(str).tolist()
        ),
        None,
    )
    if header_idx is None:
        raise ValueError("Could not locate the sampleID header row.")

    columns = raw.iloc[header_idx].tolist()
    if len(columns) != len(set(columns)):
        raise ValueError("Duplicate column names in source table.")

    df = raw.iloc[header_idx + 1:].copy()
    df.columns = columns
    df = df.dropna(subset=["sampleID"]).copy()

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Required columns missing: {sorted(missing)}")

    return df


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python scripts/01_audit_data.py <supplementary_xlsx>"
        )

    input_path = Path(sys.argv[1])
    if not input_path.is_file():
        raise SystemExit(f"Input file not found: {input_path}")

    df = load_mutations(input_path)

    for column in ["pos", "reads", "cov", "vaf", "qval"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    ref = df["ref_nt"].astype(str).str.upper()
    alt = df["mut_nt"].astype(str).str.upper()
    bases = set("ACGT")

    is_snv = (
        ref.str.len().eq(1)
        & alt.str.len().eq(1)
        & ref.isin(bases)
        & alt.isin(bases)
        & ref.ne(alt)
    )

    vaf_calculated = df["reads"] / df["cov"].where(df["cov"] != 0)
    vaf_difference = (vaf_calculated - df["vaf"]).abs()

    duplicate_count = int(
        df.duplicated(
            subset=["sampleID", "chr", "pos", "ref_nt", "mut_nt"]
        ).sum()
    )
    same_allele_count = int(ref.eq(alt).sum())

    patient_ids = df["sampleID"].astype(str).str[:7]
    patient_counts = patient_ids.value_counts().sort_index()

    report = f"""# Preliminary data audit

## Source table

- Input file: `{input_path.name}`
- Parsed records: {len(df)}
- Columns: {len(df.columns)}
- Unique sample IDs: {df['sampleID'].nunique()}
- Inferred donors (first 7 characters of sampleID): {patient_ids.nunique()}

## Data integrity checks

- Duplicate variant keys: {duplicate_count}
- Records with `ref_nt == mut_nt`: {same_allele_count}
- Missing `reads`: {int(df['reads'].isna().sum())}
- Missing or zero `cov`: {int((df['cov'].isna() | df['cov'].eq(0)).sum())}
- Missing `vaf`: {int(df['vaf'].isna().sum())}
- Rows where `vaf` differs from `reads/cov` by more than 0.00001: {int((vaf_difference > 0.00001).sum())}
- Missing `qval`: {int(df['qval'].isna().sum())}
- Records with `qval >= 0.05`: {int((df['qval'] >= 0.05).sum())}
- Missing `gene_name`: {int(df['gene_name'].isna().sum())}
- Missing `impact`: {int(df['impact'].isna().sum())}

## Variant classes

- Single-nucleotide substitutions (SNVs): {int(is_snv.sum())}
- Other records: {int((~is_snv).sum())}
- SNVs with three-character `trinuc` context: {int((is_snv & df['trinuc'].astype(str).str.len().eq(3)).sum())}

## Records per inferred donor

{patient_counts.to_string()}

## Interpretation and limitations

- The table lists variants called by ShearwaterML; it is not a raw-read dataset.
- The source documentation states that a q-value threshold of 0.05 was used.
- Records with identical reference and mutant alleles require investigation; they are not silently removed.
- Donors inferred from sample IDs are not independent biological replicates.
- This audit does not reproduce the paper's dN/dS analysis or estimate positive selection.
"""

    output_path = Path("results/data_audit.md")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report, encoding="utf-8")

    print(f"Wrote {output_path}")
    print(report)


if __name__ == "__main__":
    main()
