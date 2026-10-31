# Candidate mutation input audit

## Input summary

- Source records: 3760
- Unique biopsy IDs: 234
- Inferred donor IDs: 4
- Single-nucleotide substitutions selected: 3473
- Biopsy-level mutation observations retained: 3473

## Recurrent observations

- Donor-specific genomic events observed in more than one biopsy: 261
- Observations belonging to these recurrent events: 592
- Genomic events observed in more than one inferred donor: 74

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
