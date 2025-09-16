#!/usr/bin/env Rscript

#################################################################
source("Rscripts/functions_for_cleaning.R")

# Define directories and columns to keep
input_directories <- list(
  "h1n1pdm" = "data/raw/h1n1pdm/",
  "h3n2" = "data/raw/h3n2/", 
  "vic" = "data/raw/vic/",
  "h5n1" = "data/raw/h5n1/"
)

output_directories <- list(
  "h1n1pdm" = "data/cleaned/h1n1pdm/",
  "h3n2" = "data/cleaned/h3n2/", 
  "vic" = "data/cleaned/vic/", 
  "h5n1" = "data/cleaned/h5n1/"
)

columns_to_keep <- c("Isolate_Id", "PB2 Segment_Id", "PB1 Segment_Id", "PA Segment_Id", "HA Segment_Id", 
                     "NP Segment_Id", "NA Segment_Id", "MP Segment_Id", "NS Segment_Id", "HE Segment_Id", 
                     "P3 Segment_Id", "Isolate_Name", "Subtype", "Lineage", 
                     "Clade", "Passage_History", "Location", "Host", 
                     "Isolate_Submitter", "Submitting_Lab", "Submitting_Sample_Id", 
                     "Authors", "Publication", "Originating_Lab", "Originating_Sample_Id", 
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
# Process each directory with separate input and output directories
for (dir_name in names(input_directories)) {
  # Lineage name
  lineage <- dir_name
  
  # Step 1: Process metadata
  metadata_df <- process_metadata_files(input_directories[[dir_name]], columns_to_keep)
  counts_log <- log_counts(metadata_df, lineage, "Process Metadata", counts_log)
  
  # Step 2: Process FASTA metadata
  fasta_metadata_df <- process_fasta_files(input_directories[[dir_name]])
  
  # Step 3: Filter and combine
  filtered_df <- filter_and_combine_data(metadata_df, fasta_metadata_df, input_directories[[dir_name]])
  counts_log <- log_counts(filtered_df, lineage, "Filter and Combine", counts_log)
  
  # Step 4: Format and deduplicate
  format_dedup_df <- format_and_dedup_data(filtered_df, input_directories[[dir_name]])
  counts_log <- log_counts(format_dedup_df, lineage, "Format and Deduplicate", counts_log)
  
  # Step 5: Write outputs
  write_outputs(format_dedup_df, output_directories[[dir_name]], columns_to_keep)
}

write.csv(counts_log, paste0(Sys.Date(),"_lineage_counts_summary.csv"), row.names = FALSE)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
#### Testing script for individual strains/subtypes ####
#### USED FOR TROUBLE SHOOTING INDIVIDUAL STEPS ONLY ###
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

# lineage <- "h5n1"
# my_input_directory <- input_directories$h3n2
# my_output_directory <- output_directories$h3n2
# 
#   metadata_df <- process_metadata_files(my_input_directory, columns_to_keep)
#   
#   counts_log <- log_counts(metadata_df, lineage, "Process Metadata", counts_log)
#   
#   fasta_metadata_df <- process_fasta_files(my_input_directory)
#   
#   filtered_df <- filter_and_combine_data(metadata_df, fasta_metadata_df, my_input_directory)
#   counts_log <- log_counts(filtered_df, lineage, "Filter (missing dates) and Combine", counts_log)
#   
#   #df <- filtered_df 
#   format_dedup_df <- format_and_dedup_data(filtered_df,my_input_directory)
#   counts_log <- log_counts(format_dedup_df, lineage, "Fix Location, Passage History, Strain Name, and Deduplicate", counts_log)
#   
#   write_outputs(format_dedup_df, my_output_directory, columns_to_keep)
  