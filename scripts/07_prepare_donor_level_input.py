#!/usr/bin/env python3
"""Create donor-level SNV input for exploratory dNdScv analysis."""

from pathlib import Path
import pandas as pd

source = Path("results/dndscv_candidate_mutations.csv")
out = Path("results")
out.mkdir(exist_ok=True)

df = pd.read_csv(source)

required = {"sampleID", "donorID", "chr", "pos", "ref", "alt"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {sorted(missing)}")

event_keys = ["donorID", "chr", "pos", "ref", "alt"]

# One observation per genomic event per inferred donor.
donor_level = df.drop_duplicates(subset=event_keys).copy()

# dNdScv expects sampleID to identify the sample/donor.
donor_level = donor_level[
    ["donorID", "chr", "pos", "ref", "alt"]
].rename(columns={"donorID": "sampleID"})

donor_level.to_csv(
    out / "dndscv_donor_level_input.csv", index=False
)

# Events independently observed in more than one inferred donor.
cross_donor = (
    donor_level.groupby(["chr", "pos", "ref", "alt"])["sampleID"]
    .nunique()
)
cross_donor_count = int((cross_donor > 1).sum())

report = f"""# Donor-level dNdScv input audit

- Biopsy-level candidate observations: {len(df)}
- Inferred donors: {df["donorID"].nunique()}
- Donor-level unique mutation observations: {len(donor_level)}
- Repeated observations removed within donors: {len(df) - len(donor_level)}
- Genomic events observed in multiple inferred donors: {cross_donor_count}

Each genomic event is retained once per inferred donor. Events observed
in multiple donors are retained because they may represent independent
recurrences. Donor IDs are inferred from biopsy sample IDs and require
validation against study metadata.

This is an exploratory SNV-only input derived from a targeted dataset,
not a genome-wide mutation catalogue. Its suitability for dNdScv and
biological interpretation is therefore limited.
"""

(out / "dndscv_donor_level_audit.md").write_text(
    report, encoding="utf-8"
)

print(report)
print("Saved results/dndscv_donor_level_input.csv")
