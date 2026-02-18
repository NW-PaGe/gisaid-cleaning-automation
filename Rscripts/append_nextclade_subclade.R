#!/usr/bin/env Rscript

# Append Nextclade H5N1 clade and QC assignments to merged metadata
# Author: Pauline (WA DOH MEP)
# Purpose: Joins Nextclade clade and QC results onto the merged H5N1 metadata.xlsx
#          Run after nextclade run completes on merged/h5n1/fasta/raw_sequences_ha.fasta

library(dplyr)
library(readr)
library(openxlsx)

NEXTCLADE_TSV  <- "data/cleaned/merged/h5n1/metadata/nextclade_h5n1_ha.tsv"
METADATA_XLSX  <- "data/cleaned/merged/h5n1/metadata/metadata.xlsx"

cat("================================================================\n")
cat("    Appending Nextclade H5N1 Clade and QC Assignments\n")
cat("================================================================\n\n")

# Check inputs exist
if (!file.exists(NEXTCLADE_TSV)) {
  stop(sprintf("Nextclade TSV not found: %s", NEXTCLADE_TSV))
}
if (!file.exists(METADATA_XLSX)) {
  stop(sprintf("Merged metadata not found: %s", METADATA_XLSX))
}

# Read inputs
cat("Reading merged metadata...\n")
metadata <- read.xlsx(METADATA_XLSX)
cat(sprintf("  %d records loaded\n", nrow(metadata)))

cat("Reading Nextclade results...\n")
nextclade_results <- read_tsv(NEXTCLADE_TSV, show_col_types = FALSE)
cat(sprintf("  %d Nextclade records loaded\n", nrow(nextclade_results)))

# Verify expected columns exist - print all available cols if not to aid debugging
expected_cols <- c("seqName", "clade", "qc.overallStatus")
missing_cols <- setdiff(expected_cols, colnames(nextclade_results))
if (length(missing_cols) > 0) {
  cat("WARNING: Expected columns not found in Nextclade TSV.\n")
  cat(sprintf("  Missing: %s\n", paste(missing_cols, collapse = ", ")))
  cat(sprintf("  Available columns: %s\n", paste(colnames(nextclade_results), collapse = ", ")))
  stop("Please check column names in append_nextclade_subclade.R match the TSV output.")
}

# Select and rename relevant columns
# H5N1 Nextclade output uses 'clade' for the clade designation (e.g. 2.3.4.4b)
nextclade_clean <- nextclade_results %>%
  select(seqName, clade, qc.overallStatus) %>%
  rename(
    Isolate_Name        = seqName,
    Nextclade_Clade     = clade,
    Nextclade_QC_Status = qc.overallStatus
  )

# Drop existing Nextclade columns if re-running to avoid duplicates
metadata <- metadata %>%
  select(-any_of(c("Nextclade_Clade", "Nextclade_QC_Status")))

# Join on Isolate_Name = seqName
cat("Joining Nextclade clade and QC assignments...\n")
metadata_with_nextclade <- metadata %>%
  left_join(nextclade_clean, by = "Isolate_Name")

# Report clade assignment counts
assigned   <- sum(!is.na(metadata_with_nextclade$Nextclade_Clade))
unassigned <- sum(is.na(metadata_with_nextclade$Nextclade_Clade))
cat(sprintf("  Clade assigned:   %d sequences\n", assigned))
cat(sprintf("  Clade unassigned: %d sequences\n", unassigned))

if (unassigned > 0) {
  cat("  NOTE: Unassigned sequences may lack an HA segment or failed Nextclade QC.\n")
  cat("        Check nextclade_h5n1_ha.tsv for details.\n")
}

# Report QC status breakdown
cat("  QC status breakdown:\n")
qc_counts <- metadata_with_nextclade %>%
  count(Nextclade_QC_Status) %>%
  arrange(desc(n))
for (i in 1:nrow(qc_counts)) {
  cat(sprintf("    %-10s %d sequences\n", 
              qc_counts$Nextclade_QC_Status[i], qc_counts$n[i]))
}

# Overwrite metadata.xlsx with Nextclade columns appended
cat(sprintf("\nWriting updated metadata to %s...\n", METADATA_XLSX))
write.xlsx(metadata_with_nextclade, METADATA_XLSX, overwrite = TRUE)

cat("\nNextclade append completed successfully!\n")
