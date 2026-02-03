#!/usr/bin/env Rscript

# Load the necessary library
library(parallel)
source("Rscripts/functions_for_cleaning.R")  # Source once in main environment

# Define directories and columns to keep
input_directories <- list(
  "h1n1pdm" = "data/raw/gisaid/h1n1pdm/",
  "h3n2" = "data/raw/gisaid/h3n2/", 
  "vic" = "data/raw/gisaid/vic/",
  "h5n1" = "data/raw/gisaid/h5n1/"
)

output_directories <- list(
  "h1n1pdm" = "data/cleaned/gisaid/h1n1pdm/",
  "h3n2" = "data/cleaned/gisaid/h3n2/", 
  "vic" = "data/cleaned/gisaid/vic/", 
  "h5n1" = "data/cleaned/gisaid/h5n1/"
)

columns_to_keep <- c("Isolate_Id", "PB2 Segment_Id", "PB1 Segment_Id", "PA Segment_Id", "HA Segment_Id", 
                     "NP Segment_Id", "NA Segment_Id", "MP Segment_Id", "NS Segment_Id", "HE Segment_Id", 
                     "P3 Segment_Id", "Isolate_Name", "Subtype", "Lineage", 
                     "Clade", "Passage_History", "Location", "Host", 
                     "Isolate_Submitter", "Submitting_Sample_Id", 
                     "Publication", "Originating_Sample_Id", 
                     "Collection_Date","Note", "Update_Date", "Submission_Date", 
                     "Host_Age", "Host_Age_Unit", "Host_Gender", "Patient_Status", 
                     "Zip_Code", "Outbreak", "Human_Specimen_Source", "Animal_Specimen_Source", "Animal_Health_Status", 
                     "Domestic_Status", "PMID")

# Create a dataframe to track counts
counts_log <- data.frame(
  Lineage = character(),
  Step = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Function to process each directory
process_directory <- function(dir_name, counts_log) {
  # REMOVED: source("Rscripts/functions_for_cleaning.R")  
  # Functions are already available from parent environment
  
  lineage <- dir_name
  input_dir <- input_directories[[dir_name]]
  output_dir <- output_directories[[dir_name]]
  
  # Step 1: Process metadata
  metadata_df <- process_metadata_files(input_dir, columns_to_keep)
  counts_log <- log_counts(metadata_df, lineage, "Process Metadata", counts_log)
  
  # Step 2: Process FASTA metadata
  fasta_metadata_df <- process_fasta_files(input_dir)
  
  # Step 3: Filter and combine
  filtered_df <- filter_and_combine_data(metadata_df, fasta_metadata_df, input_dir)
  counts_log <- log_counts(filtered_df, lineage, "Filter and Combine", counts_log)
  
  # Step 4: Format and deduplicate
  format_dedup_df <- format_and_dedup_data(filtered_df, input_dir)
  counts_log <- log_counts(format_dedup_df, lineage, "Format and Deduplicate", counts_log)
  
  # Add the lineage column to format_dedup_df
  format_dedup_df$Lineage <- lineage
  
  # Step 5: Write the format_dedup_df with lineage column
  write_outputs(format_dedup_df, output_dir, columns_to_keep)
  
  # Return the processed data frame along with the updated counts_log
  return(list(format_dedup_df, counts_log))
}

# Set up the parallel cluster
num_cores <- min(detectCores() - 1, 4)  # Use max 4 cores to avoid memory issues
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
clusterExport(cl, list("input_directories", "output_directories", "columns_to_keep", 
                       "process_metadata_files", "process_fasta_files", "filter_and_combine_data", 
                       "format_and_dedup_data", "write_outputs", "log_counts", "counts_log", 
                       "process_directory", "valid_chars", "read_fasta", "is_valid_sequence",
                       "parse_fasta_header", "determine_segment_type", "process_fasta_file",
                       "filter_contemporaneous_seq", "correct_strain_format", 
                       "format_passage_dataframe", "append_egg_suffix", "check_seq_length",
                       "fix_location"))

# Load libraries in each worker
clusterEvalQ(cl, {
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(Biostrings)
})

# Use parLapply with error handling
results <- tryCatch({
  parLapply(cl, names(input_directories), function(dir_name) {
    tryCatch({
      process_directory(dir_name, counts_log)
    }, error = function(e) {
      message(sprintf("Error processing %s: %s", dir_name, conditionMessage(e)))
      return(NULL)
    })
  })
}, error = function(e) {
  message(sprintf("Fatal error in parallel processing: %s", conditionMessage(e)))
  stopCluster(cl)
  stop(e)
})

# Stop the cluster after processing
stopCluster(cl)

# Remove NULL results (failed processing)
results <- results[!sapply(results, is.null)]

if (length(results) == 0) {
  stop("All parallel processes failed. Check error messages above.")
}

# Extract the processed data frames and counts_log from the results
all_lineage_data <- lapply(results, function(res) res[[1]])
counts_log_list <- lapply(results, function(res) res[[2]])

# Combine the results from all lineages into one data frame
combined_df <- do.call(rbind, all_lineage_data)

# Write the combined data frame to a .csv file with a timestamp
output_file <- paste0("data/cleaned/gisaid_catalog/", Sys.Date(), "_combined_lineage_data.csv")
write.csv(combined_df, output_file, row.names = FALSE)

# Combine all counts_log entries and write the final counts log
counts_log <- do.call(rbind, counts_log_list)
write.csv(counts_log, paste0("data/cleaned/wa_summary_counts/", Sys.Date(), "_lineage_counts_summary.csv"), row.names = FALSE)

message("Processing completed successfully!")