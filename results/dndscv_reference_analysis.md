# Reference dNdScv analysis

## Purpose
Test the dNdScv workflow using the normal-skin example dataset bundled with the package.

## Dataset
- 3,408 mutation records
- 5 input columns: sampleID, chr, pos, ref, alt
- 74 target genes

## Results
The model returned gene-level results for all 74 target genes.
The complete output is saved in dndscv_reference_gene_results.csv.

## Limitations
This is a reference workflow test, not an analysis of the project's supplementary dataset.
The run warned about mutations at contiguous sites and mutations shared across sample IDs.
These warnings require review before interpreting results biologically.
This run should not be described as an exact reproduction of the original paper.
