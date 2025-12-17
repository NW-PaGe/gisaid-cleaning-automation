## helper functions for flu_cleaning script
# functions for cleaning 
library(tidyverse)
library(readxl)
library(writexl)
library(Biostrings)

process_metadata_files <- function(input_directories, columns_to_keep) {
  # Append the 'metadata/' sublevel to each input directory
  metadata_directories <- purrr::map(input_directories, ~ file.path(.x, "metadata"))
  
  # Process files in each metadata directory
  process_xls_file <- function(file_path, columns_to_keep) {
    data <- readxl::read_excel(file_path)
    # Only keep columns that exist in the data
    existing_columns <- intersect(columns_to_keep, names(data))
    data <- data %>%
      dplyr::select(dplyr::all_of(existing_columns))
    return(data)
  }
  
  # Process metadata files for each directory
  all_processed_data <- purrr::map(metadata_directories, function(metadata_dir) {
    # List all .xls files in the metadata directory
    xls_files <- base::list.files(metadata_dir, pattern = "*_metadata\\.xls$", full.names = TRUE)
    
    # Read and process each .xls file
    processed_data_list <- purrr::map(xls_files, ~process_xls_file(.x, columns_to_keep))
    
    # Merge dataframes for this directory
    if (base::length(processed_data_list) > 0) {
      merged_df <- purrr::reduce(processed_data_list, dplyr::full_join)
    } else {
      merged_df <- NULL  # Handle cases with no .xls files
    }
    return(merged_df)
  })
  
  # Combine all processed data from different directories into one dataframe
  final_merged_df <- purrr::reduce(all_processed_data, dplyr::full_join)
  
  return(final_merged_df)
}

# Define valid characters for nucleotide sequences
valid_chars <- c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", 
                 "a", "c", "g", "t", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v", "n")

# Function to read a FASTA file and extract valid sequences
read_fasta <- function(file_path) {
  fasta_lines <- base::readLines(file_path)
  
  # Initializing variables
  headers <- c()
  sequences <- c()
  
  current_header <- ""
  current_sequence <- ""
  
  # Loop through lines that start with > 
  for (line in fasta_lines) {
    if (base::startsWith(line, ">")) {
      if (base::nchar(current_sequence) > 0 && is_valid_sequence(current_sequence)) {
        headers <- base::c(headers, current_header)
        sequences <- base::c(sequences, current_sequence)
      }
      current_header <- base::sub("^>(.*)", "\\1", line)
      current_sequence <- ""
    } else {
      current_sequence <- base::paste0(current_sequence, line)
    }
  }
  
  # Check the last sequence
  if (base::nchar(current_sequence) > 0 && is_valid_sequence(current_sequence)) {
    headers <- base::c(headers, current_header)
    sequences <- base::c(sequences, current_sequence)
  }
  
  return(list(headers = headers, sequences = sequences))
}

# Function to check if a sequence contains only valid characters
is_valid_sequence <- function(sequence) {
  base::all(base::unlist(base::strsplit(sequence, "")) %in% valid_chars)
}

# Function to parse the header and extract Isolate_Id
parse_fasta_header <- function(header) {
  parts <- stringr::str_split(header, "\\|")[[1]]
  tibble::tibble(Isolate_Id = parts[1])
}

# Function to determine segment type based on file name
determine_segment_type <- function(file_name) {
  NA_segment <- base::grepl("_na_", file_name)
  HA_segment <- base::grepl("_ha_", file_name)
  segment <- if (NA_segment) "NA_segment" else if (HA_segment) "HA_segment" else NA
  return(list(NA_segment = NA_segment, HA_segment = HA_segment, segment = segment))
}

# Function to process a single FASTA file and return its metadata
process_fasta_file <- function(file_path) {
  fasta_data <- read_fasta(file_path)
  
  if (base::length(fasta_data$sequences) == 0) {
    return(tibble::tibble(File = base::basename(file_path), Header_ID = character(), 
                          Sequence = character(), Sequence_Length = integer(), 
                          NA_segment = logical(), HA_segment = logical(), 
                          Segment = character()))
  }
  
  metadata_list <- purrr::map_dfr(fasta_data$headers, parse_fasta_header)
  segment_info <- determine_segment_type(base::basename(file_path))
  
  metadata_list <- metadata_list %>%
    dplyr::mutate(
      Sequence = fasta_data$sequences,
      Sequence_Length = base::nchar(Sequence),  # Get sequence length
      File = base::basename(file_path),
      NA_segment = segment_info$NA_segment,
      HA_segment = segment_info$HA_segment,
      Segment = segment_info$segment
    )
  
  return(metadata_list)
}

# Main function to process all FASTA files and return combined metadata
process_fasta_files <- function(input_dir) {
  # Append 'fasta/' sublevel to the input directory
  fasta_dir <- base::file.path(input_dir, "fasta")
  
  # List all .fasta files in the fasta directory
  fasta_files <- base::list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  # Process all FASTA files and combine into one data frame
  fasta_metadata <- purrr::map_dfr(fasta_files, process_fasta_file)
  
  # Remove duplicate sequences
  fasta_metadata_dedup <- fasta_metadata %>% dplyr::distinct()
  
  return(fasta_metadata_dedup)
}

filter_and_combine_data <- function(metadata_df, fasta_metadata_df, input_directory) {
  
  # Step 1: Check for vaccine_isolate_ids.csv and load it if it exists
  vaccine_file_path <- file.path(input_directory, "vaccine_isolate_ids.csv")
  
  if (base::file.exists(vaccine_file_path)) {
    always_keep_df <- utils::read.csv(vaccine_file_path)
    isolate_ids_to_keep <- dplyr::pull(always_keep_df, Isolate_Id)
  } else {
    isolate_ids_to_keep <- base::character(0) # Empty vector if the file does not exist
  }
  
  # Step 2: Filter metadata_df with prioritization logic
  metadata_df_distinct <- metadata_df %>%
    dplyr::group_by(Isolate_Id) %>%
    dplyr::mutate(is_priority = dplyr::if_else(Isolate_Id %in% isolate_ids_to_keep, 1, 0)) %>%
    dplyr::arrange(dplyr::desc(is_priority), dplyr::desc(Update_Date), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Step 3: Ensure fasta_metadata_df is distinct
  fasta_metadata_df_distinct <- fasta_metadata_df %>%
    dplyr::select(-File) %>%
    dplyr::distinct()
  
  # Step 4: Combine metadata and fasta metadata
  complete_data_table <- metadata_df_distinct %>%
    dplyr::inner_join(fasta_metadata_df_distinct, by = "Isolate_Id")
  
  # Step 5: Filter and correct dates in the combined data
  df_filtered <- complete_data_table %>%
    dplyr::filter(!base::is.na(Collection_Date) & Collection_Date != "") %>%
    dplyr::mutate(
      Collection_Date = dplyr::case_when(
        stringr::str_detect(Collection_Date, "^\\d{4}$") ~ base::paste0(Collection_Date, "-XX-XX"),
        stringr::str_detect(Collection_Date, "^\\d{4}-\\d{2}$") ~ base::paste0(Collection_Date, "-XX"),
        TRUE ~ Collection_Date
      ),
      ha_present = dplyr::if_else(!base::is.na(`HA Segment_Id`), 1, 0),
      na_present = dplyr::if_else(!base::is.na(`NA Segment_Id`), 1, 0)
    ) %>%
    dplyr::group_by(Isolate_Id, Segment) %>%
    dplyr::mutate(is_priority = dplyr::if_else(Isolate_Id %in% isolate_ids_to_keep, 1, 0)) %>%
    dplyr::arrange(dplyr::desc(is_priority), dplyr::desc(Submission_Date), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  return(df_filtered)
}

filter_contemporaneous_seq <- function(df, input_directory) {
  # Step 1: Identify and conditionally load the vaccine_isolate_ids.csv file
  vaccine_file_path <- file.path(input_directory, "vaccine_isolate_ids.csv")
  
  if (base::file.exists(vaccine_file_path)) {
    always_keep_df <- utils::read.csv(vaccine_file_path)
    isolate_ids_to_keep <- dplyr::pull(always_keep_df, Isolate_Id)
  } else {
    isolate_ids_to_keep <- base::character(0) # Use an empty vector if the file does not exist
  }
  
  # Define passage priority: higher number = higher priority
  passage_priority <- base::c("unpassaged" = 3, "cell" = 2, "egg" = 1, "undetermined" = 0)
  
  # Step 2: Filter the dataframe
  df_filtered <- df %>%
    dplyr::mutate(
      vaccine_priority = dplyr::if_else(Isolate_Id %in% isolate_ids_to_keep, 1, 0),
      passage_priority = passage_priority[passage_category]
    ) %>%
    dplyr::arrange(Isolate_Name, Sequence, dplyr::desc(vaccine_priority), 
                   dplyr::desc(passage_priority), dplyr::desc(Submission_Date)) %>%
    dplyr::group_by(Isolate_Name, Sequence) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-vaccine_priority, -passage_priority)
  
  # Step 3: Handle contemporary isolates
  keep_contemporary_isolate_name <- df_filtered %>%
    dplyr::group_by(Isolate_Name, Segment) %>%
    dplyr::mutate(
      is_priority = dplyr::if_else(Isolate_Id %in% isolate_ids_to_keep, 1, 0),
      passage_priority = passage_priority[passage_category]
    ) %>%
    dplyr::arrange(dplyr::desc(is_priority), dplyr::desc(passage_priority), 
                   dplyr::desc(Submission_Date)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-is_priority, -passage_priority)
  
  return(keep_contemporary_isolate_name)
}


correct_strain_format <- function(df) {
  # Function to check if a strain is valid
  is_valid_strain <- function(strain) {
    if (base::grepl("QuadrivalentVaccine", strain)) return(TRUE)
    if (!base::grepl("/", strain)) return(FALSE)
    if (strain == "UnknownPassage") return(FALSE)
    
    # Regular expressions for valid strain formats
    pattern1 <- "^[A|B]/[A-Za-z-]+/([A-Za-z0-9_-]+/)*[0-9]{4}$"
    pattern2 <- "^[A|B]/[A-Za-z-]+/([A-Za-z0-9_-]+/){2}[0-9]{4}$"
    pattern3 <- "^[A|B]/([A-Za-z-]+/){2}([0-9]+/){3}[0-9]{4}$"
    
    # Check against patterns
    base::grepl(pattern1, strain) || base::grepl(pattern2, strain) || base::grepl(pattern3, strain)
  }
  
  # Function to clean Isolate_Name
  clean_isolate_name <- function(isolate_name, subtype) {
    # Standardize prefix: If it starts with "A-" or "A-/", replace with "A/"
    isolate_name <- base::gsub("^A[-/]", "A/", isolate_name)
    
    # If Subtype is B and Isolate_Name starts with A/, change it to B/
    if (subtype == "B" && base::startsWith(isolate_name, "A/")) {
      isolate_name <- base::sub("^A/", "B/", isolate_name)
    }
    
    # Remove > and ? characters
    isolate_name <- base::gsub("[>?]", "", isolate_name)
    
    # Replace multiple slashes with a single slash
    isolate_name <- base::gsub("//+", "/", isolate_name)
    
    # Change A. or A_ at the start to A/
    isolate_name <- base::gsub("^A[._]", "A/", isolate_name)
    
    # Remove all spaces in the Isolate_Name
    isolate_name <- base::gsub("\\s+", "", isolate_name)
    
    # Replace apostrophes with underscores
    isolate_name <- base::gsub("'", "_", isolate_name)
    
    # Remove diacritical marks using base R
    isolate_name <- iconv(isolate_name, from = "UTF-8", to = "ASCII//TRANSLIT")
    
    return(isolate_name)
  }
  
  # Apply cleaning to Isolate_Name
  df$Isolate_Name <- mapply(clean_isolate_name, df$Isolate_Name, df$Subtype)
  
  # Determine if strains are valid
  df$Valid_Strain <- base::sapply(df$Isolate_Name, is_valid_strain)
  
  # Separate valid and invalid strains
  valid_strains <- df[df$Valid_Strain, ]
  invalid_strains <- df[!df$Valid_Strain, ]
  
  # Combine valid and invalid strains
  final_df <- base::rbind(valid_strains, invalid_strains)
  
  # Drop the Valid_Strain column
  final_df$Valid_Strain <- NULL
  
  return(final_df)
}





# Define the function to format passage history
format_passage_dataframe <- function(df, passage_history_field = "Passage_History") {
  df <- dplyr::mutate(
    df,
    passage_category = dplyr::case_when(
      base::grepl("AM[1-9]|E[1-9]|AMNIOTIC|EGG|EX|AM_[1-9]", base::toupper(!!rlang::sym(passage_history_field))) ~ "egg",
      base::grepl("LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT", base::toupper(!!rlang::sym(passage_history_field))) ~ "unpassaged",
      base::grepl("TMK|RMK|PMK|S[1-9]|SIAT|MDCK", base::toupper(!!rlang::sym(passage_history_field))) ~ "cell",
      base::grepl("UNKNOWN|UNDEFINED", base::toupper(!!rlang::sym(passage_history_field))) ~ "undetermined",
      base::is.na(!!rlang::sym(passage_history_field)) | !!rlang::sym(passage_history_field) == "" ~ "undetermined",
      TRUE ~ "undetermined"
    )
  )
  return(df)
}

# Define the function to append an egg suffix
append_egg_suffix <- function(data) {
  data <- dplyr::mutate(
    data,
    Isolate_Name = ifelse(passage_category == "egg", paste0(Isolate_Name, "-egg"), Isolate_Name)
  )
  return(data)
}

# Define the function to check sequence length
check_seq_length <- function(df) {
  filtered_df <- dplyr::filter(
    df,
    (Segment == "HA_segment" & Sequence_Length >= 0.4 * 1700 & Sequence_Length <= 2 * 1700) |
      (Segment == "NA_segment" & Sequence_Length >= 0.4 * 1400 & Sequence_Length <= 2 * 1400)
  )
  return(filtered_df)
}

# Define the function to fix location fields
fix_location <- function(df) {
  df <- df %>%
    dplyr::mutate(
      British_Columbia_in_isolate = base::grepl("BritishColumbia", Isolate_Name),
      Idaho_in_isolate = base::grepl("Idaho", Isolate_Name),
      Washington_in_isolate = base::grepl("Washington", Isolate_Name),
      Oregon_in_isolate = base::grepl("Oregon", Isolate_Name)
    ) %>%
    dplyr::mutate(
      Location = dplyr::case_when(
        British_Columbia_in_isolate & !base::grepl("British Columbia", Location) ~ "North America / Canada / British Columbia",
        Idaho_in_isolate & !base::grepl("Idaho", Location) ~ "North America / United States / Idaho",
        Washington_in_isolate & !base::grepl("Washington", Location) ~ "North America / United States / Washington",
        Oregon_in_isolate & !base::grepl("Oregon", Location) ~ "North America / United States / Oregon",
        TRUE ~ Location
      )
    ) %>%
    dplyr::select(-British_Columbia_in_isolate, -Idaho_in_isolate, -Washington_in_isolate, -Oregon_in_isolate)
  return(df)
}

# df <- filtered_df
# input_directory <- my_input_directory

format_and_dedup_data <- function(df, input_directory){
  ## Order of functions 
  # format passage information 
  # correct strain format 
  # append egg suffix to strain name 
  # filter contemporaneous sequences 
  passage_cleaning <- format_passage_dataframe(df)
  
  correct_strain_formatting <- correct_strain_format(passage_cleaning)
  
  appending_egg <- append_egg_suffix(correct_strain_formatting)
  
  filtering_contemporaneous_sequences <- filter_contemporaneous_seq(appending_egg, input_directory)
  
  target_seq_length <- check_seq_length(filtering_contemporaneous_sequences)
  
  df_final <- fix_location(target_seq_length)
  return(df_final)
}

write_outputs <- function(df_filtered, output_dir, columns_to_keep) {
  # Filter sequences for NA and HA segments
  NA_seqs <- df_filtered %>% filter(NA_segment == TRUE)
  HA_seqs <- df_filtered %>% filter(HA_segment == TRUE)
  
  # Create DNAStringSet objects for sequences
  NA_fasta <- DNAStringSet(NA_seqs$Sequence)
  names(NA_fasta) <- NA_seqs$Isolate_Name
  
  HA_fasta <- DNAStringSet(HA_seqs$Sequence)
  names(HA_fasta) <- HA_seqs$Isolate_Name
  
  # Append 'fasta/' subdirectory to output_dir
  fasta_dir <- file.path(output_dir, "fasta")
  
  # Ensure 'fasta/' subdirectory exists, creating it if necessary
  if (!dir.exists(fasta_dir)) {
    dir.create(fasta_dir, recursive = TRUE)
  }
  
  # Write FASTA files for NA and HA sequences
  writeXStringSet(NA_fasta, filepath = file.path(fasta_dir, "raw_sequences_na.fasta"))
  writeXStringSet(HA_fasta, filepath = file.path(fasta_dir, "raw_sequences_ha.fasta"))
  
  # Append 'metadata/' subdirectory to output_dir
  metadata_dir <- file.path(output_dir, "metadata")
  
  # Ensure 'metadata/' subdirectory exists, creating it if necessary
  if (!dir.exists(metadata_dir)) {
    dir.create(metadata_dir, recursive = TRUE)
  }
  
  # Select specified columns and add additional information for metadata
  metadata <- df_filtered %>%
    select(all_of(columns_to_keep), passage_category, Sequence_Length) %>%
    distinct(Isolate_Name, .keep_all = TRUE)
  
  # Write metadata to an Excel file
  writexl::write_xlsx(metadata, file.path(metadata_dir, "metadata.xlsx"))
}

log_counts <- function(df, lineage, step, counts_log) {
  # Count the unique Isolate_Names for a specific location pattern
  count <- n_distinct(df$Isolate_Name[grepl("North America / United States / Washington", df$Location)])
  
  # Append the count to the log
  counts_log <- rbind(counts_log, data.frame(Lineage = lineage, Step = step, Count = count))
  
  return(counts_log)
}

