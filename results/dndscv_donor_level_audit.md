# Donor-level dNdScv input audit

- Biopsy-level candidate observations: 3473
- Inferred donors: 4
- Donor-level unique mutation observations: 3142
- Repeated observations removed within donors: 331
- Genomic events observed in multiple inferred donors: 74

Each genomic event is retained once per inferred donor. Events observed
in multiple donors are retained because they may represent independent
recurrences. Donor IDs are inferred from biopsy sample IDs and require
validation against study metadata.

This is an exploratory SNV-only input derived from a targeted dataset,
not a genome-wide mutation catalogue. Its suitability for dNdScv and
biological interpretation is therefore limited.
