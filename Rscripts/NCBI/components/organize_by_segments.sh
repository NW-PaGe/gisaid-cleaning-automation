#!/bin/bash

# Segment Organization Component
# Organizes sequences by segment and creates comprehensive metadata
# UPDATED: Fixes FASTA headers to use Isolate names for compatibility with GISAID

SEQUENCES_FILE="$1"    # Input FASTA file
OUTPUT_DIR="$2"        # Output directory (e.g., ncbi/h3n2)
SORT_RESULTS="$3"      # NextClade sort results (optional)

echo "    Organizing sequences by segment in $OUTPUT_DIR"

# Look for the JSON file in temp_processing
JSON_FILE=""
if [ -f "temp_processing/ncbi_dataset/data/data_report.jsonl" ]; then
    JSON_FILE="temp_processing/ncbi_dataset/data/data_report.jsonl"
elif [ -f "temp_processing/influenza_b/ncbi_dataset/data/data_report.jsonl" ]; then
    JSON_FILE="temp_processing/influenza_b/ncbi_dataset/data/data_report.jsonl"
fi

# Extract rich metadata using the new script
if [ -f "components/extract_rich_metadata.sh" ] && [ -n "$JSON_FILE" ]; then
    echo "    Using rich metadata extraction..."
    bash components/extract_rich_metadata.sh "$SEQUENCES_FILE" "$OUTPUT_DIR" "$JSON_FILE"
else
    echo "    Using simple metadata extraction (rich metadata script not found)..."
    # Fall back to simple metadata
    echo -e "accession\tfull_header\tsegment" > ${OUTPUT_DIR}/metadata/metadata.tsv
    
    grep "^>" "$SEQUENCES_FILE" | \
    while read header; do
        accession=$(echo "$header" | sed 's/^>//' | cut -d' ' -f1)
        segment=$(echo "$header" | grep -o "segment [0-9]" | grep -o "[0-9]" || echo "unknown")
        echo -e "$accession\t$header\t$segment"
    done >> ${OUTPUT_DIR}/metadata/metadata.tsv
    
    echo "      Created simple metadata for $(tail -n +2 ${OUTPUT_DIR}/metadata/metadata.tsv | wc -l) sequences"
fi

# Copy sort results if available
if [ -f "$SORT_RESULTS" ]; then
    cp "$SORT_RESULTS" ${OUTPUT_DIR}/metadata/nextclade_sort_results.tsv
    echo "      Copied NextClade sort results"
fi

# NEW: Create accession-to-isolate mapping from rich metadata
echo "    Creating accession-to-isolate name mapping..."
if [ -f "${OUTPUT_DIR}/metadata/rich_metadata.csv" ]; then
    # Extract accession and isolate name from rich metadata CSV
    # NOTE: Must use proper CSV parsing because fields like Submitters
    # contain commas (e.g. "Smith, Jones, Yang") which break naive awk parsing
    python3 -c "
import csv, sys
with open(sys.argv[1], 'r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        accession = row.get('Accession', '').strip()
        isolate = row.get('Isolate', '').strip()
        if accession and isolate:
            print(accession + '\t' + isolate)
" "${OUTPUT_DIR}/metadata/rich_metadata.csv" > /tmp/accession_to_isolate.tsv
    
    MAPPING_COUNT=$(wc -l < /tmp/accession_to_isolate.tsv)
    echo "      Created mapping for $MAPPING_COUNT accessions"
else
    echo "      WARNING: rich_metadata.csv not found, using accessions as headers"
    # Create identity mapping (accession -> accession)
    tail -n +2 ${OUTPUT_DIR}/metadata/metadata.tsv | \
    awk -F'\t' '{print $1 "\t" $1}' > /tmp/accession_to_isolate.tsv
fi

# Define segment mapping
declare -A segment_map
segment_map[1]="pb2"
segment_map[2]="pb1"
segment_map[3]="pa"
segment_map[4]="ha"
segment_map[5]="np"
segment_map[6]="na"
segment_map[7]="mp"
segment_map[8]="ns"

# Extract sequences by segment WITH HEADER RENAMING
for segment_num in {1..8}; do
    segment_name=${segment_map[$segment_num]}
    
    # Get accessions for this segment from the simple metadata (always available)
    awk -F'\t' -v seg="$segment_num" '$3==seg {print $1}' ${OUTPUT_DIR}/metadata/metadata.tsv > /tmp/${segment_name}_accessions.txt
    
    # Extract sequences
    if [ -s "/tmp/${segment_name}_accessions.txt" ]; then
        # Extract with seqtk (keeps original headers)
        seqtk subseq "$SEQUENCES_FILE" /tmp/${segment_name}_accessions.txt > /tmp/${segment_name}_temp.fasta
        
        # Rename headers to use Isolate names
        python3 << PYTHON_EOF > ${OUTPUT_DIR}/sequences/${segment_name}_sequences.fasta
import sys
import re

# Read accession-to-isolate mapping
mapping = {}
with open('/tmp/accession_to_isolate.tsv', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            accession, isolate = parts
            mapping[accession] = isolate

# Process FASTA file
with open('/tmp/${segment_name}_temp.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Extract accession (first token in header)
            header = line[1:].strip()
            accession = header.split()[0]

            # Look up isolate name from mapping
            if accession in mapping:
                isolate_name = mapping[accession]
                print(f">{isolate_name}")
            else:
                # Fallback: extract isolate from header directly
                # Handles: (A/.../2025(H3)), (A/.../2025(H3N2)), (A/.../2025)

                m = re.search(r'\(((?:A|B)/[^)]*\))\)', header)
                if m:
                    print(f">{m.group(1)}")
                else:
                    m = re.search(r'\(((?:A|B)/[^)]+)\)', header)
                    if m:
                        print(f">{m.group(1)}")
                    else:
                        # Last resort: use accession
                        print(f">{accession}")
        else:
            # Sequence line
            print(line.rstrip())
PYTHON_EOF
        
        count=$(grep -c "^>" ${OUTPUT_DIR}/sequences/${segment_name}_sequences.fasta 2>/dev/null || echo 0)
        echo "        ${segment_name}: ${count} sequences (headers renamed to Isolate names)"
        
        # Cleanup temp file
        rm -f /tmp/${segment_name}_temp.fasta
    else
        touch ${OUTPUT_DIR}/sequences/${segment_name}_sequences.fasta
        echo "        ${segment_name}: 0 sequences"
    fi
    
    # Cleanup temp file
    rm -f /tmp/${segment_name}_accessions.txt
done

# Cleanup mapping file
rm -f /tmp/accession_to_isolate.tsv

echo "      ✅ All FASTA headers renamed to match Isolate names in metadata"
