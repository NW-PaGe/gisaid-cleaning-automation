#!/usr/bin/env Rscript

# Harmonize GISAID and NCBI Data WITH DEDUPLICATION
# Author: Pauline (WA DOH MEP)
# Purpose: Merge cleaned GISAID and NCBI metadata, deduplicate by Isolate_Name

library(tidyverse)
library(readxl)
library(writexl)

cat("================================================================\n")
cat("    GISAID-NCBI Data Harmonization Pipeline\n")
cat("    WITH DEDUPLICATION BY ISOLATE NAME\n")
cat("================================================================\n\n")

# Configuration
GISAID_BASE <- "data/cleaned/gisaid"
NCBI_BASE <- "Rscripts/NCBI/ncbi"
OUTPUT_BASE <- "data/cleaned/merged"

# Create output directories
dir.create(OUTPUT_BASE, recursive = TRUE, showWarnings = FALSE)

# Define subtype mappings
subtype_mapping <- list(
  h1n1pdm = list(gisaid = "h1n1pdm", ncbi = "h1n1pdm"),
  h3n2 = list(gisaid = "h3n2", ncbi = "h3n2"),
  vic = list(gisaid = "vic", ncbi = "victoria")
)

segment_mapping <- c(
  pb2 = "PB2", pb1 = "PB1", pa = "PA", ha = "HA",
  np = "NP", na = "NA", mp = "MP", ns = "NS"
)

# Function to read GISAID metadata for a subtype
read_gisaid_metadata <- function(subtype_gisaid) {
  metadata_path <- file.path(GISAID_BASE, subtype_gisaid, "metadata", "metadata.xlsx")
  
  if (!file.exists(metadata_path)) {
    cat(sprintf("  WARNING: GISAID metadata not found for %s: %s\n", 
                subtype_gisaid, metadata_path))
    return(NULL)
  }
  
  df <- read_excel(metadata_path)
  
  # Convert ALL columns to character to ensure consistency with NCBI
  # This prevents type mismatch errors during bind_rows()
  df <- df %>%
    mutate(across(everything(), as.character))
  
  df$data_source <- "GISAID"
  return(df)
}

# Function to transform NCBI metadata to GISAID format
transform_ncbi_to_gisaid_format <- function(ncbi_df, subtype) {
  cat(sprintf("    Transforming %d NCBI segment records to GISAID format...\n", nrow(ncbi_df)))
  
  # Step 1: Create wide format - one row per isolate with segment columns
  # Group by isolate and pivot segments into columns
  ncbi_wide <- ncbi_df %>%
    mutate(
      # Convert dates to character to match GISAID format
      Collection_Date = as.character(Collection_Date),
      Release_Date = as.character(Release_Date),
      
      # Create segment column names
      Segment_Column = case_when(
        Segment == "1" ~ "PB2 Segment_Id",
        Segment == "2" ~ "PB1 Segment_Id",
        Segment == "3" ~ "PA Segment_Id",
        Segment == "4" ~ "HA Segment_Id",
        Segment == "5" ~ "NP Segment_Id",
        Segment == "6" ~ "NA Segment_Id",
        Segment == "7" ~ "MP Segment_Id",
        Segment == "8" ~ "NS Segment_Id",
        TRUE ~ NA_character_
      )
    ) %>%
    # Keep one row per isolate-segment combination
    select(
      Isolate, Accession, Segment, Segment_Column,
      Organism_Name, Submitters, Organization, Org_location,
      Release_Date, Species, Genus, Family, Host, 
      Tissue_Specimen_Source, Collection_Date, Genotype,
      Publications, Geo_Location, Country, USA
    ) %>%
    # Pivot to wide format
    pivot_wider(
      id_cols = c(Isolate, Organism_Name, Submitters, Organization, Org_location,
                  Release_Date, Species, Genus, Family, Host, 
                  Tissue_Specimen_Source, Collection_Date, Genotype,
                  Publications, Geo_Location, Country, USA),
      names_from = Segment_Column,
      values_from = Accession,
      values_fn = first  # If multiple segments, keep first
    )
  
  cat(sprintf("    Consolidated to %d unique isolates\n", nrow(ncbi_wide)))
  
  # Step 2: Transform to match GISAID column structure
  transformed <- ncbi_wide %>%
    mutate(
      # Convert ALL fields to character to match GISAID format
      # Core identifiers
      # Use the first available segment ID as Isolate_Id
      Isolate_Id = as.character(coalesce(`HA Segment_Id`, `NA Segment_Id`, `PB2 Segment_Id`, 
                            `PB1 Segment_Id`, `PA Segment_Id`, `NP Segment_Id`,
                            `MP Segment_Id`, `NS Segment_Id`)),
      Isolate_Name = as.character(Isolate),
      
      
      # Subtype and Lineage
      Subtype = as.character(case_when(
        subtype == "h3n2" ~ "A / H3N2",
        subtype == "h1n1pdm" ~ "A / H1N1",
        subtype == "vic" ~ "B / Victoria",
        TRUE ~ Genotype
      )),
      Lineage = as.character(subtype),
      
      # Location - NCBI has structured location, convert to GISAID format
      Location = as.character(case_when(
        !is.na(USA) & USA != "" ~ paste0("North America / United States / ", USA),
        Country == "USA" ~ "North America / United States",
        !is.na(Geo_Location) ~ Geo_Location,
        TRUE ~ NA_character_
      )),
      
      # Host
      Host = as.character(Host),
      
      # Submitter information
      Isolate_Submitter = as.character(Submitters),
      
      # Dates (already converted to character above)
      Update_Date = as.character(Release_Date),
      Submission_Date = as.character(Release_Date),
      
      # Additional fields (all as character to match GISAID)
      Clade = as.character(NA),
      Passage_History = as.character(NA),
      Submitting_Sample_Id = as.character(NA),
      Publication = as.character(Publications),
      Originating_Sample_Id = as.character(NA),
      Note = as.character(NA),
      Host_Age = as.character(NA),
      Host_Age_Unit = as.character(NA),
      Host_Gender = as.character(NA),
      Patient_Status = as.character(NA),
      Zip_Code = as.character(NA),
      Outbreak = as.character(NA),
      Human_Specimen_Source = as.character(Tissue_Specimen_Source),
      Animal_Specimen_Source = as.character(ifelse(Host != "Human", Tissue_Specimen_Source, NA)),
      Animal_Health_Status = as.character(NA),
      Domestic_Status = as.character(NA),
      PMID = as.character(NA),
      
      # Add data source identifier
      data_source = "NCBI"
    ) %>%
    select(
      Isolate_Id, `PB2 Segment_Id`, `PB1 Segment_Id`, `PA Segment_Id`, 
      `HA Segment_Id`, `NP Segment_Id`, `NA Segment_Id`, `MP Segment_Id`, 
      `NS Segment_Id`, Isolate_Name, Subtype, Lineage, Clade, Passage_History,
      Location, Host, Isolate_Submitter, Submitting_Sample_Id, Publication,
      Originating_Sample_Id, Collection_Date, Note, Update_Date, 
      Submission_Date, Host_Age, Host_Age_Unit, Host_Gender, Patient_Status,
      Zip_Code, Outbreak, Human_Specimen_Source, Animal_Specimen_Source,
      Animal_Health_Status, Domestic_Status, PMID, data_source
    )
  
  return(transformed)
}

# Function to read NCBI metadata for a subtype
read_ncbi_metadata <- function(subtype_ncbi) {
  metadata_path <- file.path(NCBI_BASE, subtype_ncbi, "metadata", "rich_metadata.csv")
  
  if (!file.exists(metadata_path)) {
    cat(sprintf("  WARNING: NCBI metadata not found for %s: %s\n", 
                subtype_ncbi, metadata_path))
    return(NULL)
  }
  
  df <- read_csv(metadata_path, show_col_types = FALSE)
  return(df)
}

# NEW: Function to deduplicate by Isolate_Name with priority rules
deduplicate_by_isolate_name <- function(merged_df, subtype) {
  cat(sprintf("  Deduplicating by Isolate_Name...\n"))
  
  initial_count <- nrow(merged_df)
  
  # Find duplicates
  duplicates <- merged_df %>%
    group_by(Isolate_Name) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(duplicates) > 0) {
    dup_count <- length(unique(duplicates$Isolate_Name))
    cat(sprintf("    Found %d duplicate isolate names (%d total records)\n", 
                dup_count, nrow(duplicates)))
    
    # Show examples
    cat(sprintf("    Example duplicates:\n"))
    dup_examples <- duplicates %>%
      group_by(Isolate_Name) %>%
      summarise(
        sources = paste(unique(data_source), collapse = " & "),
        count = n(),
        .groups = "drop"
      ) %>%
      head(5)
    
    for (i in 1:nrow(dup_examples)) {
      cat(sprintf("      - %s: %d records from %s\n", 
                  dup_examples$Isolate_Name[i],
                  dup_examples$count[i],
                  dup_examples$sources[i]))
    }
  } else {
    cat(sprintf("    No duplicate isolate names found\n"))
  }
  
  # Deduplication strategy:
  # 1. Prefer GISAID over NCBI (more complete metadata typically)
  # 2. Within same source, keep most recent Update_Date
  # 3. If dates equal, keep first occurrence
  
  deduplicated <- merged_df %>%
    mutate(
      # Create priority score
      source_priority = case_when(
        data_source == "GISAID" ~ 1,
        data_source == "NCBI" ~ 2,
        TRUE ~ 3
      ),
      # Parse dates for comparison (handle NA)
      update_date_parsed = as.Date(Update_Date, format = "%Y-%m-%d")
    ) %>%
    # Group by Isolate_Name
    group_by(Isolate_Name) %>%
    # Sort by priority, then date (most recent first)
    arrange(source_priority, desc(update_date_parsed), .by_group = TRUE) %>%
    # Keep first (highest priority) record
    slice(1) %>%
    ungroup() %>%
    # Remove temporary columns
    select(-source_priority, -update_date_parsed)
  
  final_count <- nrow(deduplicated)
  removed_count <- initial_count - final_count
  
  cat(sprintf("    Removed %d duplicate records\n", removed_count))
  cat(sprintf("    Final unique isolates: %d\n", final_count))
  
  # Create deduplication report
  if (removed_count > 0) {
    dedup_report <- merged_df %>%
      group_by(Isolate_Name) %>%
      filter(n() > 1) %>%
      arrange(Isolate_Name, data_source) %>%
      select(Isolate_Name, Isolate_Id, data_source, Update_Date, Collection_Date) %>%
      ungroup()
    
    report_file <- file.path(OUTPUT_BASE, subtype, "metadata", "deduplication_report.csv")
    write_csv(dedup_report, report_file)
    cat(sprintf("    Deduplication report saved to: %s\n", report_file))
  }
  
  return(deduplicated)
}

# Function to merge FASTA files with deduplication (OPTIMIZED)
merge_fasta_files <- function(gisaid_fasta, ncbi_fasta, output_fasta, deduplicated_names) {
  
  # Convert to hash set for O(1) lookup instead of O(n)
  deduplicated_set <- as.character(deduplicated_names)
  
  # Case 1: Both files exist - need to merge and deduplicate
  if (file.exists(gisaid_fasta) && file.exists(ncbi_fasta)) {
    cat(sprintf("      Merging %s...\n", basename(output_fasta)))
    
    # Use system commands for faster processing
    # Step 1: Copy GISAID file as base
    file.copy(gisaid_fasta, output_fasta, overwrite = TRUE)
    
    # Step 2: Append NCBI sequences that aren't duplicates
    # Read NCBI file and filter
    ncbi_lines <- readLines(ncbi_fasta)
    
    # Find existing GISAID headers for quick lookup
    gisaid_lines <- readLines(gisaid_fasta)
    gisaid_headers <- gisaid_lines[startsWith(gisaid_lines, ">")]
    gisaid_names <- sub("^>", "", gisaid_headers)
    
    # Parse NCBI and add non-duplicates
    to_append <- character()
    current_header <- NULL
    current_seq <- ""
    
    for (line in ncbi_lines) {
      if (startsWith(line, ">")) {
        # Save previous sequence if not duplicate
        if (!is.null(current_header)) {
          isolate_name <- sub("^>", "", current_header)
          if (!(isolate_name %in% gisaid_names)) {
            to_append <- c(to_append, current_header, current_seq)
          }
        }
        current_header <- line
        current_seq <- ""
      } else {
        current_seq <- paste0(current_seq, line)
      }
    }
    
    # Handle last sequence
    if (!is.null(current_header)) {
      isolate_name <- sub("^>", "", current_header)
      if (!(isolate_name %in% gisaid_names)) {
        to_append <- c(to_append, current_header, current_seq)
      }
    }
    
    # Append to output file
    if (length(to_append) > 0) {
      write(to_append, file = output_fasta, append = TRUE)
    }
    
    # Count sequences
    final_lines <- readLines(output_fasta)
    seq_count <- sum(startsWith(final_lines, ">"))
    return(seq_count)
    
  # Case 2: Only GISAID exists - just copy
  } else if (file.exists(gisaid_fasta)) {
    file.copy(gisaid_fasta, output_fasta, overwrite = TRUE)
    lines <- readLines(output_fasta)
    return(sum(startsWith(lines, ">")))
    
  # Case 3: Only NCBI exists - just copy
  } else if (file.exists(ncbi_fasta)) {
    file.copy(ncbi_fasta, output_fasta, overwrite = TRUE)
    lines <- readLines(output_fasta)
    return(sum(startsWith(lines, ">")))
    
  # Case 4: Neither exists - create empty file
  } else {
    writeLines(character(0), output_fasta)
    return(0)
  }
}

# Function to process a single subtype
process_subtype <- function(subtype_name, mapping) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Processing %s\n", toupper(subtype_name)))
  cat(sprintf("========================================\n"))
  
  gisaid_name <- mapping$gisaid
  ncbi_name <- mapping$ncbi
  
  # Read GISAID metadata
  cat(sprintf("  Reading GISAID %s metadata...\n", gisaid_name))
  gisaid_metadata <- read_gisaid_metadata(gisaid_name)
  
  # Read NCBI metadata
  cat(sprintf("  Reading NCBI %s metadata...\n", ncbi_name))
  ncbi_metadata <- read_ncbi_metadata(ncbi_name)
  
  if (is.null(ncbi_metadata) && is.null(gisaid_metadata)) {
    cat(sprintf("  Skipping %s - no data from either source\n", subtype_name))
    return(NULL)
  }
  
  # Transform NCBI to GISAID format if present
  if (!is.null(ncbi_metadata)) {
    cat(sprintf("  Transforming NCBI data to GISAID format...\n"))
    ncbi_transformed <- transform_ncbi_to_gisaid_format(ncbi_metadata, subtype_name)
  } else {
    ncbi_transformed <- NULL
  }
  
  # Merge metadata
  cat(sprintf("  Merging metadata...\n"))
  if (!is.null(gisaid_metadata) && !is.null(ncbi_transformed)) {
    merged_metadata <- bind_rows(gisaid_metadata, ncbi_transformed)
    cat(sprintf("    GISAID records: %d\n", nrow(gisaid_metadata)))
    cat(sprintf("    NCBI records: %d\n", nrow(ncbi_transformed)))
    cat(sprintf("    Combined (before dedup): %d\n", nrow(merged_metadata)))
  } else if (!is.null(gisaid_metadata)) {
    merged_metadata <- gisaid_metadata
    cat(sprintf("    GISAID records only: %d\n", nrow(gisaid_metadata)))
  } else {
    merged_metadata <- ncbi_transformed
    cat(sprintf("    NCBI records only: %d\n", nrow(ncbi_transformed)))
  }
  
  # DEDUPLICATE BY ISOLATE NAME
  deduplicated_metadata <- deduplicate_by_isolate_name(merged_metadata, subtype_name)
  
  # Get list of deduplicated isolate names for FASTA filtering
  deduplicated_names <- deduplicated_metadata$Isolate_Name
  
  # Create output directory
  output_dir <- file.path(OUTPUT_BASE, subtype_name)
  dir.create(file.path(output_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "fasta"), recursive = TRUE, showWarnings = FALSE)
  
  # Write deduplicated metadata
  write_xlsx(deduplicated_metadata, file.path(output_dir, "metadata", "metadata.xlsx"))
  cat(sprintf("    Written: %s\n", file.path(output_dir, "metadata", "metadata.xlsx")))
  
  # Merge FASTA files by segment WITH DEDUPLICATION
  cat(sprintf("  Merging FASTA files by segment...\n"))
  for (segment in names(segment_mapping)) {
    gisaid_fasta <- file.path(GISAID_BASE, gisaid_name, "fasta", 
                             sprintf("raw_sequences_%s.fasta", segment))
    ncbi_fasta <- file.path(NCBI_BASE, ncbi_name, "sequences", 
                           sprintf("%s_sequences.fasta", segment))
    output_fasta <- file.path(output_dir, "fasta", 
                             sprintf("raw_sequences_%s.fasta", segment))
    
    # Show which files exist for this segment
    gisaid_exists <- file.exists(gisaid_fasta)
    ncbi_exists <- file.exists(ncbi_fasta)
    
    if (!gisaid_exists && !ncbi_exists) {
      cat(sprintf("    %s: No data from either source, skipping\n", segment))
      next
    }
    
    cat(sprintf("    %s: Processing", segment))
    if (gisaid_exists && ncbi_exists) {
      cat(" (GISAID + NCBI)...")
    } else if (gisaid_exists) {
      cat(" (GISAID only)...")
    } else {
      cat(" (NCBI only)...")
    }
    
    seq_count <- merge_fasta_files(gisaid_fasta, ncbi_fasta, output_fasta, deduplicated_names)
    cat(sprintf(" %d sequences\n", seq_count))
  }
  
  return(deduplicated_metadata)
}

# Main processing
cat("\nStarting harmonization process...\n")

# Process each subtype
results <- list()
for (subtype_name in names(subtype_mapping)) {
  results[[subtype_name]] <- process_subtype(subtype_name, subtype_mapping[[subtype_name]])
}

# Generate summary report
cat("\n================================================================\n")
cat("Harmonization Complete - Summary\n")
cat("================================================================\n\n")

for (subtype_name in names(subtype_mapping)) {
  if (!is.null(results[[subtype_name]])) {
    cat(sprintf("%s:\n", toupper(subtype_name)))
    cat(sprintf("  Total sequences: %d\n", nrow(results[[subtype_name]])))
    cat(sprintf("  GISAID sequences: %d\n", 
                sum(results[[subtype_name]]$data_source == "GISAID", na.rm = TRUE)))
    cat(sprintf("  NCBI sequences: %d\n", 
                sum(results[[subtype_name]]$data_source == "NCBI", na.rm = TRUE)))
    cat(sprintf("  Output directory: %s\n\n", file.path(OUTPUT_BASE, subtype_name)))
  }
}

cat("\nNext steps:\n")
cat("1. Review deduplication reports in each subtype's metadata folder\n")
cat("2. Verify sequence counts match expectations\n")
cat("3. Use merged FASTA files for Nextstrain builds\n")
cat("4. Run phylogenetic analysis on combined datasets\n\n")

cat("Harmonization pipeline completed successfully! 🎉\n")
