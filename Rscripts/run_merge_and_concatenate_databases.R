#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readr)
library(readxl)
library(openxlsx)
library(Biostrings)

# Define paths
metadata_dir1 <- "data/cleaned/gisaid/h5n1/metadata/metadata.xlsx"
metadata_dir2 <- "data/cleaned/ncbi/h5n1/metadata/ncbi_cleaned_h5n1_wide.csv"
merged_metadata_dir <- "data/cleaned/merged/h5n1/metadata/"
fasta_dir1 <- "data/cleaned/gisaid/h5n1/fasta/raw_sequences_ha.fasta"
fasta_dir2 <- "data/cleaned/ncbi/h5n1/fasta/raw_sequences_ha.fasta"
merged_fasta_dir <- "data/cleaned/merged/h5n1/fasta/"

# Create output directories if they do not exist
dir.create(merged_metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_fasta_dir, recursive = TRUE, showWarnings = FALSE)

# Step 1: Read and merge metadata files
cat("Reading metadata files...\n")
metadata1 <- read_excel(metadata_dir1) %>% mutate(Source = "GISAID")
metadata2 <- read_csv(metadata_dir2) %>% mutate(Source = "NCBI")

# Ensure Collection_Date is a Date
#metadata1$Collection_Date <- as.Date(metadata1$Collection_Date, format = "%Y-%m-%d")
#metadata2$Collection_Date <- as.Date(metadata2$Collection_Date, format = "%Y-%m-%d")
# Collection_Date needs to be treated as a character to allow for the XX months and days

# Function to extract FASTA headers and sequences
tidy_fasta <- function(fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  data.frame(
    Isolate_Name = names(fasta),
    Sequence = as.character(fasta),
    stringsAsFactors = FALSE
  )
}

# Extract FASTA data
cat("Processing FASTA files...\n")
fasta_data1 <- tidy_fasta(fasta_dir1)
fasta_data2 <- tidy_fasta(fasta_dir2)

# Merge FASTA data with corresponding metadata
cat("Merging FASTA data with metadata...\n")
metadata1 <- left_join(metadata1, fasta_data1, by = "Isolate_Name")
metadata2 <- left_join(metadata2, fasta_data2, by = "Isolate_Name")

# Merge metadata1 and metadata2, keeping NCBI sequences when duplicate Isolate_Names exist
cat("Merging metadata files and resolving duplicates...\n")
merged_metadata <- bind_rows(metadata1, metadata2)

merged_metadata_final <- merged_metadata %>%
  mutate(Source = factor(Source, levels = c("NCBI", "GISAID"))) %>%
  arrange(desc(Source)) %>%  # If found in GISAID then take GISAID entry
  distinct(Isolate_Name, .keep_all = TRUE)  # Keeps first occurrence per Isolate_Name

# Write merged metadata to XLSX
output_metadata_file <- file.path(merged_metadata_dir, "metadata.xlsx")
write.xlsx(merged_metadata_final, output_metadata_file, overwrite = TRUE)
cat("Merged metadata saved to:", output_metadata_file, "\n")

# Merge FASTA sequences, keeping only the selected (deduplicated) metadata entries
filtered_fasta <- merged_metadata_final %>%
  filter(!is.na(Sequence)) %>%
  select(Isolate_Name, Sequence)

# Convert to DNAStringSet
total_fasta <- DNAStringSet(filtered_fasta$Sequence)
names(total_fasta) <- filtered_fasta$Isolate_Name

# Write merged FASTA file
output_fasta_file <- file.path(merged_fasta_dir, "raw_sequences_ha.fasta")
writeXStringSet(total_fasta, output_fasta_file)
cat("Merged FASTA saved to:", output_fasta_file, "\n")

cat("Processing completed successfully!\n")
