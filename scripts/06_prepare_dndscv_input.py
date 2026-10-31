#!/usr/bin/env python3
"""Prepare and audit candidate SNVs while preserving biopsy-level observations."""

from pathlib import Path
import sys

import pandas as pd


def load_table(path):
    raw = pd.read_excel(path, header=None)

    header_idx = next(
        i for i in range(len(raw))
        if "sampleID" in raw.iloc[i].astype(str).tolist()
    )

    df = raw.iloc[header_idx + 1:].copy()
    df.columns = raw.iloc[header_idx].tolist()

    return df.dropna(subset=["sampleID"]).copy()


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python scripts/06_prepare_dndscv_input.py <input.xlsx>"
        )

    source = Path(sys.argv[1])
    df = load_table(source)

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["ref_nt"] = df["ref_nt"].astype(str).str.upper()
    df["mut_nt"] = df["mut_nt"].astype(str).str.upper()
    df["donorID"] = df["sampleID"].astype(str).str[:7]

    bases = set("ACGT")

    snv = (
        df["ref_nt"].str.len().eq(1)
        & df["mut_nt"].str.len().eq(1)
        & df["ref_nt"].isin(bases)
        & df["mut_nt"].isin(bases)
        & df["ref_nt"].ne(df["mut_nt"])
        & df["pos"].notna()
        & df["chr"].notna()
    )

    candidate = df.loc[
        snv, ["sampleID", "donorID", "chr", "pos", "ref_nt", "mut_nt"]
    ].copy()

    candidate.columns = [
        "sampleID", "donorID", "chr", "pos", "ref", "alt"
    ]

    event_keys = ["chr", "pos", "ref", "alt"]

    # Count the number of biopsies in which each event was observed
    # within each inferred donor.
    recurrence = (
        candidate.groupby(["donorID"] + event_keys)["sampleID"]
        .nunique()
        .reset_index(name="n_biopsies")
    )

    recurrent = recurrence.loc[
        recurrence["n_biopsies"] > 1
    ].copy()

    # Preserve all biopsy-level observations in the candidate table.
    # Do not collapse repeated events at this stage.
    out = Path("results")
    out.mkdir(parents=True, exist_ok=True)

    candidate.to_csv(
        out / "dndscv_candidate_mutations.csv", index=False
    )

    recurrent.to_csv(
        out / "recurrent_mutations_within_donors.csv", index=False
    )

    repeated_event_count = len(recurrent)
    recurrent_observation_count = int(
        candidate.merge(
            recurrent[["donorID"] + event_keys],
            on=["donorID"] + event_keys,
            how="inner"
        ).shape[0]
    )

    cross_donor = (
        candidate[["donorID"] + event_keys]
        .drop_duplicates()
        .groupby(event_keys)["donorID"]
        .nunique()
    )

    cross_donor_event_count = int((cross_donor > 1).sum())

    report = f"""# Candidate mutation input audit

## Input summary

- Source records: {len(df)}
- Unique biopsy IDs: {df["sampleID"].nunique()}
- Inferred donor IDs: {df["donorID"].nunique()}
- Single-nucleotide substitutions selected: {int(snv.sum())}
- Biopsy-level mutation observations retained: {len(candidate)}

## Recurrent observations

- Donor-specific genomic events observed in more than one biopsy: {repeated_event_count}
- Observations belonging to these recurrent events: {recurrent_observation_count}
- Genomic events observed in more than one inferred donor: {cross_donor_event_count}

## Interpretation and limitations

The original biopsy-level sample IDs are retained. Repeated observations
within a donor are not automatically treated as independent mutations
or removed, because they may reflect spatially distributed clonal
expansions.

Donor IDs are inferred from the first seven characters of sample IDs
and should be checked against the study metadata.

This table is a candidate mutation dataset, not a validated dNdScv input.
The suitability of the targeted sequencing design, sample structure,
and mutation calls for dNdScv must be assessed before model fitting
or biological interpretation.
"""

    (out / "dndscv_input_audit.md").write_text(
        report, encoding="utf-8"
    )

    print(report)
    print(f"Wrote {out / 'dndscv_candidate_mutations.csv'}")
    print(f"Wrote {out / 'recurrent_mutations_within_donors.csv'}")


if __name__ == "__main__":
    main()
