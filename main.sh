#!/bin/bash
set -e  # Exit on any error

echo "Starting data cleaning process..."

# Get S3 bucket paths from environment variables (with fallbacks for safety)
S3_INPUT_BUCKET="${S3_INPUT_BUCKET:-s3://your-input-bucket/data/raw/}"
S3_OUTPUT_BUCKET="${S3_OUTPUT_BUCKET:-s3://your-output-bucket/data/cleaned/}"
LOCAL_RAW_DIR="data/raw/"
LOCAL_CLEANED_DIR="data/cleaned/"

# R scripts to run in order
RSCRIPT1="Rscripts/run_main_cleaning_parallel.R"
RSCRIPT2="Rscripts/run_clean_ncbi_data_parallel.R"
RSCRIPT3="Rscripts/run_merge_and_concatenate_databases.R"

# Navigate to the root directory
cd "$(dirname "$0")"

# Display which buckets we're using (for logging/debugging)
echo "Input bucket: $S3_INPUT_BUCKET"
echo "Output bucket: $S3_OUTPUT_BUCKET"

# Sync new data from S3
echo "Syncing new data from S3..."
aws s3 sync "$S3_INPUT_BUCKET" "$LOCAL_RAW_DIR"

# Remove the cleaned data directory
echo "Removing local cleaned data folder..."
rm -rf "$LOCAL_CLEANED_DIR"

# Create required directory structure
echo "Creating required directory structure..."
BASE_DIR="data/cleaned"
REQUIRED_DIRS=(
    "gisaid/h1n1pdm/fasta"
    "gisaid/h1n1pdm/metadata"
    "gisaid/h3n2/fasta"
    "gisaid/h3n2/metadata"
    "gisaid/vic/fasta"
    "gisaid/vic/metadata"
    "gisaid/h5n1/fasta"
    "gisaid/h5n1/metadata"
    "gisaid_catalog"
    "merged/h5n1/fasta"
    "merged/h5n1/metadata"
    "ncbi/h5n1/fasta"
    "ncbi/h5n1/metadata"
    "wa_summary_counts"
)

for DIR in "${REQUIRED_DIRS[@]}"; do
    mkdir -p "${BASE_DIR}/${DIR}"
done

echo "All required directories have been created."

# Run R scripts in sequence
echo "Starting R script execution..."
for RSCRIPT in "$RSCRIPT1" "$RSCRIPT2" "$RSCRIPT3"; do
    if [ -f "$RSCRIPT" ]; then
        echo "Running R script: $RSCRIPT"
        Rscript "$RSCRIPT"
        echo "Completed: $RSCRIPT"
    else
        echo "Error: R script '$RSCRIPT' not found. Exiting."
        exit 1
    fi
done

echo "All R scripts executed successfully!"

# Upload cleaned data back to S3
echo "Uploading cleaned data to S3..."
aws s3 sync "$LOCAL_CLEANED_DIR" "$S3_OUTPUT_BUCKET"

echo "Data cleaning process completed successfully!"
