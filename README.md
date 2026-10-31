# Somatic Mutation Analysis in Normal Human Skin

**A reproducible bioinformatics analysis of somatic variation, mutational patterns and preliminary signals of positive selection in normal human epidermis.**

## Project overview

Somatic mutations accumulate in normal tissues throughout life. In the skin, ultraviolet exposure and tissue-specific selection contribute to the development and expansion of somatically mutated cell populations. Characterising these mutations provides insight into the mutational processes operating in healthy tissues and the early stages of clonal evolution.

This project investigates somatic mutation calls reported in normal human skin by Martincorena et al. (2015). It implements a reproducible Python workflow for data quality assessment, variant classification, mutational spectrum analysis and gene-level annotation summaries. It also includes preparation of donor-level mutation data and a preliminary dN/dS analysis using **dNdScv**.

The analysis was designed to address three related questions:

1. What are the principal characteristics and quality-control properties of the available somatic variant dataset?
2. What substitution patterns and gene-level mutation distributions are observed in the reported calls?
3. Can the available mutation data provide an informative preliminary assessment of selection using a dN/dS framework?

The descriptive analyses were completed, and dNdScv was successfully executed on a candidate set of 3,142 donor-level SNVs. The resulting estimates suggest an excess of nonsynonymous substitutions under the fitted model. However, the analysis remains exploratory because the source dataset is targeted, the independence of recurrent events requires further assessment, and the suitability of the input for dNdScv requires additional validation.

The project therefore distinguishes observed mutation patterns from model-based evidence of selection.

---

## Data source

The analysis uses Supplementary Dataset S1 from:

Martincorena I, Roshan A, Gerstung M, et al. (2015). **High burden and pervasive positive selection of somatic mutations in normal human skin.** *Science*, 348(6237), 880–886.

- [Publication DOI](https://doi.org/10.1126/science.aaa6806)
- [Full text – PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC4471149/)

The supplementary workbook contains somatic variant calls generated using ShearwaterML. Available annotations include:

- Genomic coordinates and reference alleles
- Alternate alleles
- Read counts and sequencing coverage
- Variant allele fractions
- Trinucleotide sequence contexts
- Gene annotations
- Predicted functional impact
- Protein and nucleotide change annotations

The source table contains 3,760 records from 234 biopsy sample IDs, representing four inferred donors.

**Data availability and scope**

The source is a targeted somatic mutation dataset, not a genome-wide catalogue of all mutations and not a raw sequencing dataset. Consequently, the analyses are restricted to the reported calls and annotations. They do not involve read-level variant calling or independent validation of the original mutation calls.

The supplementary workbook is not redistributed with this repository. To reproduce the analysis, obtain the supplementary material from the original publication and place `NIHMS63718-supplement-2.xlsx` in the `data/` directory.

---

## Repository structure

```text
normal-skin-positive-selection-starter/
│
├── data/
│   └── NIHMS63718-supplement-2.xlsx
│
├── scripts/
│   ├── 01_audit_data.py
│   ├── 02_mutational_spectrum.py
│   ├── 03_classify_variants.py
│   ├── 04_check_uv_dinucleotide.py
│   ├── 05_gene_impact_summary.py
│   └── 06_prepare_dndscv_input.py
│
├── results/
│   ├── data_audit.md
│   ├── variant_classification.md
│   ├── mutational_spectrum.md
│   ├── uv_dinucleotide_check.md
│   ├── gene_impact_summary.md
│   ├── dndscv_input_audit.md
│   ├── dndscv_reference_analysis.md
│   ├── dndscv_candidate_mutations.csv
│   ├── dndscv_donor_level_input.csv
│   ├── own_dndscv_global_results.csv
│   ├── own_dndscv_gene_results.csv
│   ├── own_dndscv_full_results.rds
│   └── [additional CSV tables and figures]
│
├── environment.yml
└── README.md
```

---

## Computational environment

### Python

The Python workflow was developed using Python 3.11 and the following principal packages:

- pandas
- openpyxl
- matplotlib
- seaborn

Create and activate the Conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate normal-skin-positive-selection
```

### R

The dNdScv analysis was performed separately using:

- R 4.4.3
- dndscv 0.0.0.9
- IRanges

The `IRanges` package was explicitly loaded before running dNdScv.

The Python and R analyses use separate environments.

---

## Reproducing the analysis

Run the Python scripts from the repository root:

```bash
python scripts/01_audit_data.py data/NIHMS63718-supplement-2.xlsx

python scripts/02_mutational_spectrum.py data/NIHMS63718-supplement-2.xlsx

python scripts/03_classify_variants.py data/NIHMS63718-supplement-2.xlsx

python scripts/04_check_uv_dinucleotide.py data/NIHMS63718-supplement-2.xlsx

python scripts/05_gene_impact_summary.py data/NIHMS63718-supplement-2.xlsx

python scripts/06_prepare_dndscv_input.py data/NIHMS63718-supplement-2.xlsx
```

The scripts generate quality-control reports, summary tables and figures in `results/`.

The dNdScv analysis uses the donor-level candidate table generated during the workflow:

```text
results/dndscv_donor_level_input.csv
```

Example R command used for the analysis:

```r
library(IRanges)
library(dndscv)

mutations <- read.csv(
  "results/dndscv_donor_level_input.csv",
  stringsAsFactors = FALSE
)

results <- dndscv(mutations)

write.csv(
  results$globaldnds,
  "results/own_dndscv_global_results.csv",
  row.names = FALSE
)

write.csv(
  results$sel_cv,
  "results/own_dndscv_gene_results.csv",
  row.names = FALSE
)

saveRDS(
  results,
  "results/own_dndscv_full_results.rds"
)
```

---

## Analysis and results

### 1. Data quality assessment

The initial audit examined record counts, sample identifiers, duplicate variant keys, missing annotations, sequencing coverage, read counts and variant allele fractions.

The supplementary table contains:

| Metric | Result |
|---|---:|
| Total variant records | 3,760 |
| Unique biopsy sample IDs | 234 |
| Inferred donors | 4 |
| Duplicate records under the audit variant key | 0 |
| Records with identical reference and alternate alleles | 2 |
| Missing gene annotations | 251 |
| Missing impact annotations | 257 |

No missing read counts, coverage values or variant allele fractions were identified.

Reported variant allele fractions were consistent with the ratio of alternate reads to total coverage (`reads / cov`), within a tolerance of 0.00001 for all records.

Donor identifiers were inferred from the first seven characters of the sample identifiers. This inference was used to organise the downstream analysis and should be validated against the original study's sampling metadata.

**Output:** `results/data_audit.md`

### 2. Variant classification

Variants were classified according to the lengths and contents of their reference and alternate allele fields.

| Variant class | Records |
|---|---:|
| Single-nucleotide substitutions (SNVs) | 3,473 |
| Equal-length multi-nucleotide changes | 203 |
| Deletion-like representations | 64 |
| Unequal-length sequence changes | 18 |
| Identical reference and alternate alleles | 2 |
| **Total** | **3,760** |

This classification describes the representation of variants in the supplementary table. It does not replace reference-genome-based normalization or independent validation of individual calls.

**Outputs:**
- `results/variant_classification.md`
- `results/variant_class_counts.csv`

### 3. Single-nucleotide substitution spectrum

The substitution spectrum was calculated from 3,473 SNVs.

All selected SNVs had valid trinucleotide contexts, with the central nucleotide matching the reported reference allele.

Substitutions were represented using pyrimidine-oriented notation:

| Substitution | Count |
|---|---:|
| C>T | 1,981 |
| C>A | 536 |
| T>A | 323 |
| T>C | 314 |
| T>G | 180 |
| C>G | 139 |
| **Total** | **3,473** |

The analysis also generated a 96-channel trinucleotide-context summary, accounting for the two flanking reference nucleotides.

The observed spectrum is descriptive. No sequence-opportunity normalisation or formal mutational-signature fitting was performed. Therefore, the spectrum alone cannot establish a specific mutational signature or quantify the contribution of ultraviolet radiation.

**Outputs:**
- `results/mutational_spectrum.md`
- `results/six_class_counts.csv`
- `results/96_context_counts.csv`
- `results/mutational_spectrum_96.png`

### 4. Multi-nucleotide changes and potential UV-associated patterns

Equal-length multi-nucleotide changes were examined separately from single-nucleotide substitutions.

Among 203 equal-length multi-nucleotide records, the most frequent allele patterns were:

| Allele change | Count |
|---|---:|
| GG>AA | 98 |
| CC>TT | 59 |

These patterns are relevant to the study of mutagenesis in sun-exposed skin. However, the presence of these allele changes does not by itself establish their biological origin.

A separate inspection of adjacent SNV calls identified one adjacent pair in the same original biopsy, reported as a two-base change in `NOTCH1` (`C1264T,C1265T`).

Adjacent genomic positions in a donor-level aggregated table were not treated as confirmed multi-nucleotide events. Such proximity alone is insufficient to establish that two substitutions arose as one event.

Further interpretation would require reference-genome validation and assessment of the original variant representation.

**Outputs:**
- `results/uv_dinucleotide_check.md`
- `results/multi_nucleotide_change_counts.csv`

### 5. Gene-level variant impact summary

Gene-level summaries were generated using the available gene and predicted impact annotations.

A total of 3,503 records had both annotations available, representing 75 genes. Records missing either annotation were excluded from this summary.

The genes with the highest numbers of annotated records included:

| Gene | Annotated records |
|---|---:|
| NOTCH1 | 451 |
| NOTCH2 | 194 |
| FAT1 | 147 |
| MUC17 | 147 |
| FAT4 | 143 |
| APOB | 133 |

Impact categories included missense, synonymous, nonsense, essential splice-site, frameshift deletion and in-frame deletion variants.

These counts reflect the composition of the targeted supplementary dataset. They are not normalised for gene length, sequence context, mutation opportunity or sample-specific effects. Consequently, the number of variants in a gene should not be interpreted as a direct measure of positive selection.

**Outputs:**
- `results/gene_impact_summary.md`
- `results/gene_impact_counts.csv`

### 6. Preparation of donor-level mutation data

A candidate SNV table was prepared for downstream dN/dS analysis.

The original dataset contains multiple biopsies per inferred donor. To reduce repeated observations of the same genomic event across biopsies from a donor, SNVs were deduplicated using the following fields:

- Donor identifier
- Chromosome
- Genomic position
- Reference allele
- Alternate allele

| Processing step | Records |
|---|---:|
| Original supplementary records | 3,760 |
| Selected SNVs | 3,473 |
| Unique donor-level SNVs after deduplication | 3,142 |
| Repeated biopsy-level observations removed | 331 |
| Genomic events observed in more than one donor | 74 |
| Observations belonging to these cross-donor events | 152 |

The resulting donor-level table contains 3,142 mutation records assigned to four donor identifiers.

The deduplication procedure is intended to reduce potential pseudoreplication arising from repeated observations of mutations in multiple biopsies from the same donor. It does not establish clonal relationships or prove that every repeated observation represents the same biological event.

Events observed in different donors were retained because they may represent independent recurrent mutations. Their independence has not been experimentally established.

The resulting table is therefore a candidate input rather than a fully validated representation of independent mutational events.

**Outputs:**
- `results/dndscv_input_audit.md`
- `results/dndscv_candidate_mutations.csv`
- `results/dndscv_donor_level_input.csv`
- `results/recurrent_mutations_within_donors.csv`

### 7. Preliminary dN/dS analysis using dNdScv

The candidate donor-level dataset was analysed using dNdScv to obtain preliminary estimates of the ratio of nonsynonymous to synonymous substitutions.

The analysis was performed on 3,142 donor-level SNVs across four inferred donors.

#### Global estimates

| Mutation class | dN/dS estimate | 95% confidence interval |
|---|---:|---:|
| Missense | 1.427 | 1.181–1.724 |
| Nonsense | 1.716 | 1.178–2.500 |
| Splice | 2.021 | 1.197–3.413 |
| Truncating | 1.805 | 1.308–2.492 |
| All substitutions | 1.461 | 1.213–1.760 |

The global estimate for all substitutions was 1.46 (95% CI: 1.21–1.76), indicating an excess of nonsynonymous substitutions relative to the model's expectation.

#### Gene-level results

The following genes had `qallsubs_cv` values below 0.05 in the fitted model:

| Gene | Synonymous | Missense | Nonsense | Splice | q-value (all substitutions) |
|---|---:|---:|---:|---:|---:|
| AURKA | 0 | 4 | 0 | 1 | 0.0078 |
| SUFU | 0 | 3 | 0 | 0 | 0.0094 |
| CDKN2A | 0 | 4 | 0 | 0 | 0.0114 |
| NRAS | 0 | 3 | 0 | 0 | 0.0331 |

These results identify candidate gene-level signals for further investigation. However, the mutation counts in several of these genes are low, and the corresponding dN/dS estimates are highly sensitive to small changes in the observed counts.

For example, very large missense dN/dS estimates in genes with only three or four observed missense mutations should not be interpreted as precise measures of selection strength.

The analysis also produced a warning that the same genomic mutations were present in multiple donor IDs. These observations may represent independent recurrent events, but this assumption requires further evaluation.

**Important methodological qualification**

The source dataset is targeted and does not constitute an unbiased genome-wide catalogue of somatic mutations. The candidate input also contains SNVs only and does not reproduce the complete variant composition of the original study's dNdScv example dataset.

The assumptions of dNdScv, including the suitability of the mutation catalogue, sequence context, genomic reference and independence of observations, require further assessment before the gene-level estimates can be considered definitive evidence of positive selection.

The dNdScv results are therefore presented as preliminary model-based findings rather than a final biological conclusion.

**Outputs:**
- `results/own_dndscv_global_results.csv`
- `results/own_dndscv_gene_results.csv`
- `results/own_dndscv_full_results.rds`

A separate reference workflow was also tested using the normal-skin example dataset bundled with the dNdScv package. This was a software workflow test and is distinct from the analysis of the project's candidate mutation table.

---

## Limitations

Several limitations should be considered when interpreting the results:

1. **Targeted sequencing design.** The source dataset represents mutations reported in selected genes and cannot be assumed to provide an unbiased genome-wide mutation catalogue.

2. **Donor assignment.** Donor identifiers were inferred from sample ID prefixes. This relationship should be checked against the original study metadata.

3. **Clonal structure.** Multiple biopsies from the same donor may contain related cell populations. Deduplication reduces repeated observations but does not reconstruct clonal relationships.

4. **Recurrent mutations.** Variants observed in multiple donors may represent independent events at recurrent hotspots. Their independence has not been established.

5. **Multi-nucleotide variants.** Allele patterns and adjacent SNV positions require reference-genome validation before they can be confidently classified as single multi-nucleotide events.

6. **dNdScv input suitability.** The candidate input contains SNVs selected from a targeted dataset. Its compatibility with the assumptions of the dNdScv model requires further evaluation.

7. **Low mutation counts.** Some gene-level estimates are based on few observed mutations and may be unstable.

8. **Mutational signatures.** No formal signature decomposition or sequence-opportunity-normalised analysis was performed.

These limitations restrict the strength of the biological conclusions and define the next steps required for a more rigorous selection analysis.

---

## Conclusion

This project establishes a reproducible workflow for examining somatic mutation calls in normal human skin, from initial data quality assessment to preliminary dN/dS modelling.

The descriptive analyses characterise the substitution spectrum, multi-nucleotide allele patterns and distribution of annotated variants across genes. A donor-level candidate mutation table was subsequently constructed and analysed using dNdScv.

The preliminary model returned a global dN/dS estimate above one and identified several candidate gene-level signals. These findings are consistent with an excess of nonsynonymous substitutions under the fitted model, but their biological interpretation remains conditional on validation of the input data and model assumptions.

The work provides a foundation for further investigation of somatic evolution in normal human tissue, particularly the relationship between mutational processes, recurrent driver mutations and positive selection.

---

## References

1. Martincorena I, Roshan A, Gerstung M, et al. (2015). High burden and pervasive positive selection of somatic mutations in normal human skin. *Science*, 348(6237), 880–886. https://doi.org/10.1126/science.aaa6806

2. Martincorena I, Raine KM, Gerstung M, et al. (2017). Universal Patterns of Selection in Cancer and Somatic Tissues. *Cell*, 171(5), 1029–1041.e21. https://doi.org/10.1016/j.cell.2017.09.042

3. dNdScv: dN/dS estimation and selection analysis. https://github.com/martincorena/dndscv

---

## Reproducibility and data availability

The repository contains the analysis scripts, environment specification, generated reports, summary tables and model outputs.

The original supplementary workbook is not included and must be obtained from the publication. All reported results refer to the specific supplementary dataset and processing steps described above.

The analyses are intended to be reproducible from the repository root, subject to obtaining the source workbook and installing the specified software dependencies.