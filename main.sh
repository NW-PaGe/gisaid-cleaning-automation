#!/bin/bash
set -e  # Exit on any error

echo "Starting data cleaning process..."

# Get S3 bucket paths from environment variables (with fallbacks for safety)
S3_INPUT_BUCKET="${S3_INPUT_BUCKET:-s3://your-input-bucket/data/raw/}"
S3_OUTPUT_BUCKET="${S3_OUTPUT_BUCKET:-s3://your-output-bucket/data/cleaned/}"
LOCAL_RAW_DIR="data/raw/"
LOCAL_CLEANED_DIR="data/cleaned/"

# R scripts to run in order (keep intact)
RSCRIPT1="Rscripts/run_main_cleaning_parallel.R"
RSCRIPT2="Rscripts/run_clean_ncbi_data_parallel.R"
RSCRIPT3="Rscripts/run_merge_and_concatenate_databases.R"

# NEW: NCBI seasonal pipeline + harmonization
NCBI_SEASONAL_PIPELINE="Rscripts/NCBI/waphl_main_pipeline.sh"
HARMONIZE_SCRIPT="Rscripts/harmonize_gisaid_ncbi.R"

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

# Create required directory structure (keep intact)
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
    "ncbi/h1n1pdm/sequences"
    "ncbi/h1n1pdm/metadata"
    "ncbi/h3n2/sequences"
    "ncbi/h3n2/metadata"
    "ncbi/vic/sequences"
    "ncbi/vic/metadata"
    "wa_summary_counts"
)

for DIR in "${REQUIRED_DIRS[@]}"; do
    mkdir -p "${BASE_DIR}/${DIR}"
done

echo "All required directories have been created."

# Run R scripts in sequence (KEEP INTACT)
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

# -------------------------------------------------------------------
# NEW: Run seasonal NCBI pipeline (bash) from Rscripts/NCBI
# -------------------------------------------------------------------
echo ""
echo "Starting NCBI seasonal pipeline (bash): $NCBI_SEASONAL_PIPELINE"
if [ -f "$NCBI_SEASONAL_PIPELINE" ]; then
    pushd "Rscripts/NCBI" >/dev/null
    chmod +x "./waphl_main_pipeline.sh"
    ./waphl_main_pipeline.sh
    popd >/dev/null
    echo "Completed NCBI seasonal pipeline."
else
    echo "Error: NCBI seasonal pipeline script '$NCBI_SEASONAL_PIPELINE' not found. Exiting."
    exit 1
fi

# -------------------------------------------------------------------
# Copy NCBI seasonal results into data/cleaned/ncbi/ so they are
# picked up by the harmonization script and the final S3 sync
# -------------------------------------------------------------------
echo ""
echo "Copying NCBI seasonal results to data/cleaned/ncbi/..."
for subtype in h1n1pdm h3n2 vic; do
    src="Rscripts/NCBI/ncbi/${subtype}"
    dest="${LOCAL_CLEANED_DIR}ncbi/${subtype}"
    if [ -d "$src" ]; then
        cp -r "$src"/. "$dest"/
        echo "  Copied ${subtype}: $(find "$dest" -name '*.fasta' -o -name '*.csv' -o -name '*.tsv' | wc -l) files"
    else
        echo "  WARNING: $src not found, skipping"
    fi
done
echo "NCBI seasonal results copied."

# -------------------------------------------------------------------
# NEW: Harmonize GISAID + NCBI outputs (run from repo root)
# -------------------------------------------------------------------
echo ""
echo "Starting harmonization (R): $HARMONIZE_SCRIPT"
if [ -f "$HARMONIZE_SCRIPT" ]; then
    Rscript "$HARMONIZE_SCRIPT"
    echo "Completed harmonization."
else
    echo "Error: Harmonization script '$HARMONIZE_SCRIPT' not found. Exiting."
    exit 1
fi

# Upload cleaned data back to S3 (MOVED TO END)
echo ""
echo "Uploading cleaned data to S3..."
aws s3 sync "$LOCAL_CLEANED_DIR" "$S3_OUTPUT_BUCKET"

echo "Data cleaning process completed successfully!"
