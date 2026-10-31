#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dndscv)
  library(IRanges)
  library(IRanges)
})

input <- "results/dndscv_donor_level_input.csv"
out <- "results"

if (!file.exists(input)) {
  stop("Missing donor-level input. Run script 07 first.")
}

mutations <- read.csv(input, stringsAsFactors = FALSE)

required <- c("sampleID", "chr", "pos", "ref", "alt")
if (!all(required %in% names(mutations))) {
  stop("Input file is missing required dNdScv columns.")
}

if (anyDuplicated(mutations[c("sampleID", "chr", "pos", "ref", "alt")])) {
  stop("Duplicate donor-level mutation records detected.")
}

results <- dndscv(mutations)

write.csv(
  results$globaldnds,
  file.path(out, "own_dndscv_global_results.csv"),
  row.names = FALSE
)

write.csv(
  results$sel_cv,
  file.path(out, "own_dndscv_gene_results.csv"),
  row.names = FALSE
)

saveRDS(
  results,
  file.path(out, "own_dndscv_full_results.rds")
)

cat("dNdScv analysis completed.\n")
cat("Global results and gene-level results saved in results/.\n")
print(results$globaldnds)
