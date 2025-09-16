#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
# Function to read and concatenate metadata files with renaming, transformations, and segment conversion
read_metadata_files <- function(directory) {
  metadata_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  
  metadata <- metadata_files %>%
    lapply(function(file) {
      read_csv(file, col_types = cols(
        Release_Date = col_character(),   
        Collection_Date = col_character()
      )) %>%
        mutate(
          # Trim whitespace and standardize format
          Release_Date = str_trim(Release_Date),
          Collection_Date = str_trim(Collection_Date),
          Collection_Date = str_replace_all(Collection_Date, "/", "-"),
          
          # Format Collection_Date to YYYY-MM-DD with "XX" placeholders
          Collection_Date = case_when(
            str_detect(Collection_Date, "^\\d{4}$") ~ paste0(Collection_Date, "-XX-XX"),  # "2024" → "2024-XX-XX"
            str_detect(Collection_Date, "^\\d{4}-\\d{2}$") ~ paste0(Collection_Date, "-XX"),  # "2024-05" → "2024-05-XX"
            str_detect(Collection_Date, "^\\d{4}-\\d{2}-\\d{2}$") ~ Collection_Date,  # Keep full dates unchanged
            TRUE ~ NA_character_  # Keep other unknown formats as NA
          ),
          
          # Format Release_Date the same way
          Release_Date = case_when(
            str_detect(Release_Date, "^\\d{4}$") ~ paste0(Release_Date, "-XX-XX"),
            str_detect(Release_Date, "^\\d{4}-\\d{2}$") ~ paste0(Release_Date, "-XX"),
            str_detect(Release_Date, "^\\d{4}-\\d{2}-\\d{2}$") ~ Release_Date,
            TRUE ~ NA_character_
          )
        ) 
    }) %>%
    bind_rows() %>%
    
    # Rename variables
    rename(
      Isolate_Id = Accession,
      Isolate_Submitter = Submitters,
      Submitting_Lab = Organization
    ) %>%
    
    # Create Location variable
    mutate(Location = case_when(
      Geo_Location %in% c("USA: Washington", "USA: WA") ~ "North America / United States / Washington",
      TRUE ~ Geo_Location
    )) %>%
    
    # Create Subtype and Lineage variables
    mutate(
      Subtype = if_else(Genotype == "H5N1", "A / H5N1", NA_character_),
      Lineage = if_else(Genotype == "H5N1", "h5n1", NA_character_)
    ) %>%
    
    # Extract Isolate_Name from GenBank_Title and clean up trailing subtype info
    mutate(Isolate_Name = str_extract(GenBank_Title, "\\(([^)]+)\\)")) %>%  
    mutate(Isolate_Name = str_remove_all(Isolate_Name, "[()]")) %>%        
    mutate(Isolate_Name = str_remove(Isolate_Name, "H[0-9]+N[0-9]+$")) %>%  
    
    # Convert Segment numbers to their short abbreviations
    mutate(Segment_Abbreviation = sapply(Segment, convert_segment_to_abbreviation))
  
  return(metadata)
}


# Function to convert segment numbers to abbreviations
convert_segment_to_abbreviation <- function(segment_number) {
  # Lookup table for segment numbers to abbreviations
  segment_lookup <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
  
  # Convert the segment number to its abbreviation, accounting for 1-based indexing
  if (segment_number >= 1 && segment_number <= 8) {
    return(segment_lookup[segment_number])
  } else {
    return(NA)  # Return NA if the segment number is not valid
  }
}

# Function to read and concatenate FASTA files if there are multiples
read_fasta_files <- function(directory) {
  # List all FASTA files in the directory
  fasta_files <- list.files(directory, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in the directory.")
  }
  
  # Process each FASTA file
  fasta <- lapply(fasta_files, function(file) {
    lines <- readLines(file)
    
    # Remove any carriage return characters
    lines <- gsub("\r", "", lines)
    
    # Identify headers (lines starting with ">") and group sequences
    is_header <- grepl("^>", lines)
    groups <- cumsum(is_header)  # Create groups for headers and their corresponding sequences
    
    # Combine sequences under each header
    fasta_data <- split(lines, groups) %>%
      lapply(function(block) {
        header <- block[1]  # First line is the header
        seq <- paste(block[-1], collapse = "")  # Concatenate multi-line sequences
        data.frame(Header = header, Sequence = seq, stringsAsFactors = FALSE)
      })
    
    do.call(rbind, fasta_data)
  }) %>%
    bind_rows()
  
  # Trim any extra whitespace from headers and sequences
  fasta <- fasta %>%
    mutate(across(everything(), str_trim))  # Trim whitespace
  
  # Separate the 'Header' column into 'Isolate_Id' and 'GenBank_Title_Fasta' using regex
  fasta <- fasta %>%
    mutate(Header = str_remove(Header, "^>"),  # Remove the leading '>'
           Isolate_Id = str_extract(Header, "^[^|]+"),  # Everything before the first '|'
           GenBank_Title = str_remove(Header, "^[^|]+\\|"))  # Everything after the first '|'
  
  # Add a column for sequence length
  fasta <- fasta %>%
    mutate(Sequence_Length = nchar(Sequence))  # Count the number of characters in the sequence
  
  # Check for mismatched headers and sequences
  mismatched <- fasta %>%
    filter(is.na(Sequence) | Sequence == "")
  
  if (nrow(mismatched) > 0) {
    warning(paste("Mismatch between headers and sequences in the following entries:",
                  paste(mismatched$Isolate_Id, collapse = ", ")))
  }
  
  # Select and reorder columns
  fasta <- fasta %>%
    select(Isolate_Id, GenBank_Title, Sequence, Sequence_Length)
  
  return(fasta)
}

# Function to merge metadata and FASTA files
merge_metadata_and_fasta <- function(metadata, fasta) {
  metadata <- metadata %>%
    mutate(Isolate_Id = str_trim(Isolate_Id)) %>% 
    mutate(GenBank_Title = str_trim(GenBank_Title))
  
  fasta <- fasta %>%
    mutate(Isolate_Id = str_trim(Isolate_Id)) %>% 
    mutate(GenBank_Title = str_trim(GenBank_Title))
  
  merged <- metadata %>%
    full_join(fasta, by = c("Isolate_Id","GenBank_Title"))
  return(merged)
}

# Function to export sequences by segment to separate FASTA files
export_sequences_by_segment <- function(cleaned_data, output_directory) {
  # Ensure the output directory exists
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Iterate through each segment abbreviation (e.g., HA, NA, etc.)
  unique_segments <- unique(cleaned_data$Segment_Abbreviation)
  
  for (segment in unique_segments) {
    # Filter the data for the current segment
    segment_data <- cleaned_data %>%
      filter(Segment_Abbreviation == segment)
    
    # Create a FASTA header and sequence for each isolate
    fasta_content <- segment_data %>%
      mutate(Fasta_Header = paste(">", Isolate_Name, sep = ""),
             Fasta_Sequence = Sequence) %>%
      select(Fasta_Header, Fasta_Sequence) %>%
      pmap_chr(~ paste(.x, .y, sep = "\n"))  # Combine header and sequence
    
    # Write the sequences to a FASTA file
    output_filepath <- file.path(output_directory, paste0("raw_sequences_", tolower(segment), ".fasta"))
    writeLines(fasta_content, output_filepath)
    
    # Print a message indicating progress
    message(paste("FASTA file created for segment", segment, ":", output_filepath))
  }
}
# Function to create the wide format from the metadata
create_wide_format_metadata <- function(metadata) {
  wide_metadata <- metadata %>%
    select(
      Organism_Name, Release_Date, Species, Genus, Family, Host,
      Tissue_Specimen_Source, Collection_Date, Location, Subtype, 
      Lineage, Isolate_Name, Segment_Abbreviation, Isolate_Id, Isolate_Submitter, Submitting_Lab
    ) %>%
    # Pivot the data to create a wide format, where each Segment_Abbreviation becomes a new column
    pivot_wider(
      names_from = Segment_Abbreviation,
      values_from = Isolate_Id,
      names_glue = "{Segment_Abbreviation} Segment_Id"  # Add space between segment abbreviation and '_Segment_Id'
    )
  
  return(wide_metadata)
}

# Main script function updated to include the wide format
main <- function() {
  # Define input and output directories
  metadata_dir <- "data/raw/ncbi/h5n1/metadata/"
  fasta_dir <- "data/raw/ncbi/h5n1/fasta/"
  output_dir <- "data/cleaned/ncbi/h5n1/metadata"
  output_fasta_dir <- "data/cleaned/ncbi/h5n1/fasta"  # New output directory for segment-specific FASTA files
  
  # Read metadata and FASTA files
  metadata <- read_metadata_files(metadata_dir)
  fasta <- read_fasta_files(fasta_dir)
  
  # Merge metadata and FASTA files
  cleaned_data <- merge_metadata_and_fasta(metadata, fasta)
  # Clean Isolate_Name by removing spaces and special characters
  cleaned_data <- cleaned_data %>%
    mutate(Isolate_Name = str_replace_all(Isolate_Name, "[' ]", ""))
  # Create and save the wide-format metadata
  wide_metadata <- create_wide_format_metadata(cleaned_data)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(wide_metadata, file.path(output_dir, "ncbi_cleaned_h5n1_wide.csv"), row.names = FALSE)
  
  # Save the regular cleaned data (long format)
  write.csv(cleaned_data, file.path(output_dir, "ncbi_cleaned_h5n1.csv"), row.names = FALSE)
  
  # Export sequences by segment to separate FASTA files
  export_sequences_by_segment(cleaned_data, output_fasta_dir)
}

# Run the main script
main()

