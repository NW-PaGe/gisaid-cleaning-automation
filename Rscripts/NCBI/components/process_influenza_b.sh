#!/bin/bash

# Component: Process Influenza B (Victoria lineage)
# Purpose: Download, deduplicate, and organize WAPHL Influenza B sequences with isolate-based FASTA headers

set -euo pipefail

START_DATE=${START_DATE:-"2025-01-01"}
OUTPUT_BASE=${OUTPUT_BASE:-"ncbi"}

echo ""
echo "========================================="
echo "PROCESSING INFLUENZA B"
echo "========================================="
echo ""

# Create temp directory for influenza B data
mkdir -p temp_ncbi_data_b

echo "Step 1B: Downloading all Influenza B sequences from Washington..."
# Influenza B virus taxon ID: 11520
datasets download virus genome taxon 11520 \
  --geo-location "USA: Washington" \
  --released-after "${START_DATE}" \
  --filename temp_ncbi_data_b/influenza_b_wa.zip

echo "Step 2B: Extracting Influenza B download..."
cd temp_ncbi_data_b
unzip -q -o influenza_b_wa.zip
cd ..

# Count total sequences
TOTAL_SEQUENCES_B=$(grep -c "^>" temp_ncbi_data_b/ncbi_dataset/data/genomic.fna 2>/dev/null || echo "0")
echo "Total Influenza B sequences downloaded: ${TOTAL_SEQUENCES_B}"

if [ "${TOTAL_SEQUENCES_B}" -eq 0 ]; then
  echo "WARNING: No Influenza B sequences found!"
  exit 0
fi

echo "Step 3B: Finding WAPHL Influenza B sequence accessions..."
# Extract accessions for sequences with WAPHL in headers
grep "^>" temp_ncbi_data_b/ncbi_dataset/data/genomic.fna | \
  grep -i "WAPHL" | \
  sed 's/^>//' | \
  cut -d' ' -f1 > temp_ncbi_data_b/waphl_accessions.txt

# If no WAPHL found, try broader WA DOH search
if [ ! -s temp_ncbi_data_b/waphl_accessions.txt ]; then
  echo "No WAPHL sequences found, trying broader Washington State DOH search..."
  grep "^>" temp_ncbi_data_b/ncbi_dataset/data/genomic.fna | \
    grep -i "washington.*health\|health.*washington\|WA.*PHL\|WADOH" | \
    sed 's/^>//' | \
    cut -d' ' -f1 > temp_ncbi_data_b/waphl_accessions.txt
fi

WAPHL_COUNT_B=$(wc -l < temp_ncbi_data_b/waphl_accessions.txt)
echo "Found ${WAPHL_COUNT_B} WAPHL Influenza B sequences"

if [ "${WAPHL_COUNT_B}" -eq 0 ]; then
  echo "WARNING: No WAPHL sequences identified in Influenza B data"
  exit 0
fi

echo "Step 4B: Extracting WAPHL Influenza B sequences using seqtk..."
seqtk subseq temp_ncbi_data_b/ncbi_dataset/data/genomic.fna \
  temp_ncbi_data_b/waphl_accessions.txt > temp_ncbi_data_b/waphl_sequences.fna

echo "Step 5B: Deduplicating Influenza B sequences..."
# Simple deduplication - keep first occurrence of each header line
seqtk seq temp_ncbi_data_b/waphl_sequences.fna | \
  awk '/^>/ {if (seen[$0]++) next} {print}' > temp_ncbi_data_b/waphl_sequences_dedup.fna

DEDUP_COUNT=$(grep -c "^>" temp_ncbi_data_b/waphl_sequences_dedup.fna 2>/dev/null || echo "0")
echo "  After deduplication: ${DEDUP_COUNT} sequences"

# Output directory for Victoria
mkdir -p "${OUTPUT_BASE}/vic"/{metadata,sequences}

echo "Step 6B: Extracting rich Influenza B metadata..."
if [ -f "components/extract_rich_metadata.sh" ]; then
  echo "  Using rich metadata extraction script..."
  bash components/extract_rich_metadata.sh \
    temp_ncbi_data_b/waphl_sequences_dedup.fna \
    "${OUTPUT_BASE}/vic" \
    temp_ncbi_data_b/ncbi_dataset/data/data_report.jsonl

  METADATA_COUNT=$(tail -n +2 "${OUTPUT_BASE}/vic/metadata/rich_metadata.csv" 2>/dev/null | wc -l || echo "0")
  echo "  Rich metadata records extracted: ${METADATA_COUNT}"

  if [ "${METADATA_COUNT}" -eq 0 ]; then
    echo "WARNING: No metadata records extracted for Influenza B!"
  fi
else
  echo "WARNING: components/extract_rich_metadata.sh not found; organize step may fall back to simple metadata"
fi

echo "Step 7B: Organizing Influenza B sequences by segment (with header renaming)..."
if [ ! -f "components/organize_by_segments.sh" ]; then
  echo "ERROR: components/organize_by_segments.sh not found; cannot organize by segment"
  exit 1
fi

# IMPORTANT:
# - Do NOT re-run seqtk subseq per segment after this.
# - organize_by_segments.sh already writes per-segment FASTAs AND renames headers to isolate.
bash components/organize_by_segments.sh \
  temp_ncbi_data_b/waphl_sequences_dedup.fna \
  "${OUTPUT_BASE}/vic" \
  ""

# Summarize counts from the produced per-segment FASTAs
VIC_TOTAL=0
echo "  Segment FASTA counts (post-rename):"
for f in "${OUTPUT_BASE}/vic/sequences/"*_sequences.fasta; do
  [ -e "$f" ] || continue
  c=$(grep -c "^>" "$f" 2>/dev/null || echo "0")
  echo "    $(basename "$f"): $c"
  VIC_TOTAL=$((VIC_TOTAL + c))
done

echo ""
echo "✅ Influenza B processing complete!"
echo "   Total sequences (sum across segments): ${VIC_TOTAL}"
echo "   Sequences organized in: ${OUTPUT_BASE}/vic/sequences/"
echo "   Metadata available in:  ${OUTPUT_BASE}/vic/metadata/"
echo ""