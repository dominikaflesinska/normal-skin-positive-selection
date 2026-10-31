# Preliminary data audit

## Source table

- Input file: `NIHMS63718-supplement-2.xlsx`
- Parsed records: 3760
- Columns: 19
- Unique sample IDs: 234
- Inferred donors (first 7 characters of sampleID): 4

## Data integrity checks

- Duplicate variant keys: 0
- Records with `ref_nt == mut_nt`: 2
- Missing `reads`: 0
- Missing or zero `cov`: 0
- Missing `vaf`: 0
- Rows where `vaf` differs from `reads/cov` by more than 0.00001: 0
- Missing `qval`: 0
- Records with `qval >= 0.05`: 0
- Missing `gene_name`: 251
- Missing `impact`: 257

## Variant classes

- Single-nucleotide substitutions (SNVs): 3473
- Other records: 287
- SNVs with three-character `trinuc` context: 3473

## Records per inferred donor

sampleID
PD13634     891
PD18003    1676
PD20399     962
PD21910     231

## Interpretation and limitations

- The table lists variants called by ShearwaterML; it is not a raw-read dataset.
- The source documentation states that a q-value threshold of 0.05 was used.
- Records with identical reference and mutant alleles require investigation; they are not silently removed.
- Donors inferred from sample IDs are not independent biological replicates.
- This audit does not reproduce the paper's dN/dS analysis or estimate positive selection.
