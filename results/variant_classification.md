# Variant classification

## Dataset

- Source file: `NIHMS63718-supplement-2.xlsx`
- Total records: 3760
- Number of classified records: 3760
- Classification totals match the input record count: True

## Variant classes

variant_class
Single-nucleotide substitution               3473
Equal-length multi-nucleotide change          203
Deletion-like representation                   64
Unequal-length sequence change                 18
Identical reference and alternate alleles       2

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
