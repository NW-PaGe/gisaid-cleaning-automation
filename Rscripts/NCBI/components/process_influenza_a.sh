#!/bin/bash

# Influenza A Processing Component (H3N2, H1N1pdm)
# Part of WAPHL NCBI Pipeline

set -e

START_DATE="$1"
OUTPUT_BASE="ncbi"

echo "Processing Influenza A sequences (H3N2, H1N1pdm)..."
echo "Start date: $START_DATE"

# Create directory structure
mkdir -p ${OUTPUT_BASE}/{h3n2,h1n1pdm}/{sequences,metadata}
mkdir -p temp_processing

# Install seqtk if needed
if ! command -v seqtk &> /dev/null; then
    echo "  Installing seqtk..."
    mamba install -c bioconda seqtk -y
fi

# Download Influenza A
echo "  Downloading Influenza A from NCBI..."
datasets download virus genome taxon 197911 \
  --geo-location "USA: Washington" \
  --released-after ${START_DATE} \
  --filename temp_processing/influenza_a_wa.zip

# Extract sequences
cd temp_processing
unzip -q influenza_a_wa.zip
cd ..

TOTAL_SEQUENCES=$(grep -c "^>" temp_processing/ncbi_dataset/data/genomic.fna)
echo "  Total sequences: ${TOTAL_SEQUENCES}"

# Extract WAPHL sequences
echo "  Extracting WAPHL sequences..."
grep "^>" temp_processing/ncbi_dataset/data/genomic.fna | \
grep -i "WAPHL\|washington.*health" | \
sed 's/^>//' | cut -d' ' -f1 > temp_processing/waphl_accessions.txt

WAPHL_COUNT=$(wc -l < temp_processing/waphl_accessions.txt)
echo "  WAPHL sequences found: ${WAPHL_COUNT}"

if [ "$WAPHL_COUNT" -eq 0 ]; then
    echo "  ❌ No WAPHL sequences found"
    exit 1
fi

# Extract complete sequences
seqtk subseq temp_processing/ncbi_dataset/data/genomic.fna temp_processing/waphl_accessions.txt > temp_processing/waphl_sequences.fasta

# Use NextClade sort to classify by subtype
echo "  Running NextClade sort for subtype classification..."
nextclade sort \
  --output-dir temp_processing/sorted \
  --output-results-tsv temp_processing/sort_results.tsv \
  temp_processing/waphl_sequences.fasta

# Process sorted results
echo "  Processing sorted sequences..."

# Function to process a subtype
process_subtype() {
    local subtype_lower=$1
    local subtype_upper=$2
    local search_pattern=$3
    
    echo "    Processing ${subtype_upper} sequences..."
    
    # Find directories matching the pattern
    found_dir=""
    for dir in temp_processing/sorted/nextstrain/flu/${search_pattern}* temp_processing/sorted/*${search_pattern}*; do
        if [ -d "$dir" ]; then
            found_dir="$dir"
            break
        fi
    done
    
    if [ -n "$found_dir" ]; then
        echo "      Found ${subtype_upper} directory: $found_dir"
        
        # Combine all FASTA files from this directory
        find "$found_dir" -name "*.fasta" -exec cat {} \; > temp_processing/${subtype_lower}_all.fasta
        
        if [ -s temp_processing/${subtype_lower}_all.fasta ]; then
            # Deduplicate sequences
            awk '
            /^>/ {
                split($1, parts, ".");
                accession = parts[1];
                gsub(">", "", accession);
                if (!seen[accession]) {
                    seen[accession] = 1;
                    print_seq = 1;
                    print $0;
                } else {
                    print_seq = 0;
                }
                next;
            }
            print_seq { print $0 }
            ' temp_processing/${subtype_lower}_all.fasta > temp_processing/${subtype_lower}_dedup.fasta
            
            # Organize by segments
            bash components/organize_by_segments.sh temp_processing/${subtype_lower}_dedup.fasta ${OUTPUT_BASE}/${subtype_lower} temp_processing/sort_results.tsv
            
            dedup_count=$(grep -c "^>" temp_processing/${subtype_lower}_dedup.fasta)
            echo "      ${subtype_upper} sequences after deduplication: $dedup_count"
        else
            echo "      No ${subtype_upper} sequences found"
        fi
    else
        echo "      No ${subtype_upper} directory found in NextClade sort output"
    fi
}

# Process H3N2
process_subtype "h3n2" "H3N2" "h3n2"

# Process H1N1pdm  
process_subtype "h1n1pdm" "H1N1pdm" "h1n1"

echo "✓ Influenza A processing complete"
exit 0
