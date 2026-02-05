#!/bin/bash

# Rich Metadata Extraction Component  
# Extracts comprehensive metadata from NCBI datasets JSON to match H5N1 CSV format
# FIXED: Properly extracts isolate names from FASTA headers as fallback

SEQUENCES_FILE="$1"    # Input FASTA file
OUTPUT_DIR="$2"        # Output directory (e.g., ncbi/h3n2)
JSON_FILE="$3"         # NCBI datasets JSON file (data_report.jsonl)

echo "    Extracting rich metadata to match H5N1 CSV format..."

# Create the enhanced metadata extraction script
cat > /tmp/extract_rich_metadata.py << 'EOF'
#!/usr/bin/env python3

import json
import csv
import sys
import re
from collections import defaultdict

def extract_isolate_from_fasta_header(header: str):
    """
    Correctly extracts isolate from headers like:
      Influenza A virus (A/.../2025(H3))      -> A/.../2025(H3)
      Influenza A virus (A/.../2025(H3N2))    -> A/.../2025(H3N2)
      Influenza A virus (A/.../2025)          -> A/.../2025
    Same for B/.
    """

    # 1) Double-paren case: (A/.../2025(H3)) or (A/.../2025(H3N2))
    m = re.search(r'\(((?:A|B)/[^)]*\))\)', header)
    if m:
        return m.group(1)  # includes the inner closing paren

    # 2) Single-paren case: (A/.../2025)
    m = re.search(r'\(((?:A|B)/[^)]+)\)', header)
    if m:
        return m.group(1)

    return None

def extract_rich_metadata(json_file, output_file, sequences_file):
    """Extract rich metadata from NCBI datasets JSON to match H5N1 CSV format"""
    
    # First, get isolate names from FASTA headers
    print("Reading FASTA headers to extract isolate names...", file=sys.stderr)
    fasta_isolates = {}  # accession -> isolate_name
    fasta_titles = {}    # accession -> full title
    
    try:
        with open(sequences_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    parts = header.split()
                    if parts:
                        accession = parts[0]
                        
                        # Extract isolate name
                        isolate = extract_isolate_from_fasta_header(header)
                        if isolate:
                            fasta_isolates[accession] = isolate
                        
                        # Store full title (everything after accession)
                        if len(parts) > 1:
                            fasta_titles[accession] = ' '.join(parts[1:])
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
    
    print(f"Extracted isolate names for {len(fasta_isolates)} sequences from FASTA", file=sys.stderr)
    
    # Get list of sequence accessions
    sequence_accessions = set(fasta_isolates.keys())
    print(f"Processing {len(sequence_accessions)} sequences", file=sys.stderr)
    
    # Define the CSV columns matching H5N1 format
    fieldnames = [
        'Accession', 'Organism_Name', 'GenBank_RefSeq', 'Assembly', 'SRA_Accession',
        'Submitters', 'Organization', 'Org_location', 'Release_Date', 'Isolate',
        'Species', 'Genus', 'Family', 'Molecule_type', 'Length', 'Nuc_Completeness',
        'Genotype', 'Segment', 'Publications', 'Geo_Location', 'Country', 'USA',
        'Host', 'Tissue_Specimen_Source', 'Collection_Date', 'BioSample', 'BioProject',
        'GenBank_Title'
    ]
    
    metadata_records = []
    processed_count = 0
    
    try:
        with open(json_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    data = json.loads(line.strip())
                    accession = data.get('accession', '')
                    
                    # Only process sequences that are in our FASTA file
                    if accession not in sequence_accessions:
                        continue
                    
                    processed_count += 1
                    
                    # Extract data with safe navigation
                    def safe_get(obj, *keys):
                        """Safely navigate nested dictionary"""
                        current = obj
                        for key in keys:
                            if isinstance(current, dict) and key in current:
                                current = current[key]
                            else:
                                return ''
                        return current if current is not None else ''
                    
                    def get_lineage_name(lineage_list, level):
                        """Extract specific taxonomic level from lineage"""
                        if not isinstance(lineage_list, list):
                            return ''
                        for item in lineage_list:
                            if isinstance(item, dict) and item.get('name', '').lower() == level.lower():
                                return item.get('name', '')
                        return ''
                    
                    # Extract submitter names
                    submitter_names = safe_get(data, 'submitter', 'names')
                    if isinstance(submitter_names, list):
                        submitters = ', '.join(submitter_names)
                    else:
                        submitters = ''
                    
                    # Extract organism name
                    organism_name = safe_get(data, 'virus', 'organismName')
                    
                    # Extract taxonomic info from virus lineage
                    virus_lineage = safe_get(data, 'virus', 'lineage')
                    family = get_lineage_name(virus_lineage, 'Orthomyxoviridae')
                    genus = get_lineage_name(virus_lineage, 'Alphainfluenzavirus')
                    species = safe_get(data, 'virus', 'organismName')
                    
                    # Extract host info
                    host_name = safe_get(data, 'host', 'organismName')
                    
                    # Extract geographic info
                    geo_location = safe_get(data, 'location', 'geographicLocation')
                    usa_state = safe_get(data, 'location', 'usaState')
                    country = 'USA' if usa_state else safe_get(data, 'submitter', 'country')
                    
                    # Extract dates
                    collection_date = safe_get(data, 'isolate', 'collectionDate')
                    release_date = safe_get(data, 'releaseDate')
                    if release_date:
                        # Convert ISO format to date only
                        release_date = release_date.split('T')[0]
                    
                    # CRITICAL FIX: Get isolate name from FASTA header first, fallback to JSON
                    isolate_name = fasta_isolates.get(accession, '')
                    if not isolate_name:
                        # Fallback to JSON
                        isolate_name = safe_get(data, 'isolate', 'name')
                    
                    # If still empty, try to construct from organism name
                    if not isolate_name and organism_name:
                        isolate_name = organism_name.replace('Influenza A virus ', '').replace('Influenza B virus ', '')
                    
                    # Determine genotype from isolate name (basic extraction)
                    genotype = ''
                    if isolate_name:
                        if '(H3N2)' in isolate_name or 'H3N2' in isolate_name or '/h3n2' in isolate_name.lower():
                            genotype = 'H3N2'
                        elif '(H1N1)' in isolate_name or 'H1N1' in isolate_name or '/h1n1' in isolate_name.lower():
                            genotype = 'H1N1'
                        elif '(H5N1)' in isolate_name or 'H5N1' in isolate_name or '/h5n1' in isolate_name.lower():
                            genotype = 'H5N1'
                        elif 'Victoria' in isolate_name or 'Yamagata' in isolate_name:
                            genotype = 'Victoria' if 'Victoria' in isolate_name else 'Yamagata'
                    
                    # Create the metadata record
                    record = {
                        'Accession': accession,
                        'Organism_Name': organism_name,
                        'GenBank_RefSeq': 'GenBank',  # All from GenBank in our case
                        'Assembly': '',  # Not typically available in viral data
                        'SRA_Accession': '',  # Would need separate lookup
                        'Submitters': submitters,
                        'Organization': safe_get(data, 'submitter', 'affiliation'),
                        'Org_location': safe_get(data, 'submitter', 'country'),
                        'Release_Date': release_date,
                        'Isolate': isolate_name,
                        'Species': species,
                        'Genus': genus,
                        'Family': family,
                        'Molecule_type': 'ssRNA(-)',  # Standard for influenza
                        'Length': str(safe_get(data, 'length')),
                        'Nuc_Completeness': safe_get(data, 'completeness').lower(),
                        'Genotype': genotype,
                        'Segment': str(safe_get(data, 'segment')),
                        'Publications': '',  # Would need separate lookup
                        'Geo_Location': geo_location,
                        'Country': country,
                        'USA': usa_state,
                        'Host': host_name,
                        'Tissue_Specimen_Source': safe_get(data, 'isolate', 'source'),
                        'Collection_Date': collection_date,
                        'BioSample': '',  # Would need separate lookup
                        'BioProject': '',  # Would need separate lookup
                        'GenBank_Title': fasta_titles.get(accession, '')
                    }
                    
                    metadata_records.append(record)
                    
                except json.JSONDecodeError as e:
                    print(f"Error parsing JSON line {line_num}: {e}", file=sys.stderr)
                    continue
                except Exception as e:
                    print(f"Error processing line {line_num}: {e}", file=sys.stderr)
                    continue
    
    except FileNotFoundError:
        print(f"Error: JSON file {json_file} not found", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading JSON file: {e}", file=sys.stderr)
        return
    
    # Write CSV output
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(metadata_records)
        
        print(f"Successfully processed {processed_count} metadata records", file=sys.stderr)
        
        # Show sample of isolate names
        print("Sample isolate names:", file=sys.stderr)
        for record in metadata_records[:5]:
            print(f"  {record['Accession']}: {record['Isolate']}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error writing CSV file: {e}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_rich_metadata.py <json_file> <output_csv> <sequences_fasta>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    output_csv = sys.argv[2]
    sequences_fasta = sys.argv[3]
    
    extract_rich_metadata(json_file, output_csv, sequences_fasta)
EOF

# Run the metadata extraction
if [ -f "$JSON_FILE" ] && [ -f "$SEQUENCES_FILE" ]; then
    echo "      Extracting rich metadata from JSON and FASTA..."
    python3 /tmp/extract_rich_metadata.py "$JSON_FILE" "${OUTPUT_DIR}/metadata/rich_metadata.csv" "$SEQUENCES_FILE"
    
    if [ -f "${OUTPUT_DIR}/metadata/rich_metadata.csv" ]; then
        echo "      âœ“ Rich metadata extracted to rich_metadata.csv"
        
        # Also create the simple TSV for compatibility
        echo -e "accession\tfull_header\tsegment" > ${OUTPUT_DIR}/metadata/metadata.tsv
        
        grep "^>" "$SEQUENCES_FILE" | \
        while read header; do
            accession=$(echo "$header" | sed 's/^>//' | cut -d' ' -f1)
            segment=$(echo "$header" | grep -o "segment [0-9]" | grep -o "[0-9]" || echo "unknown")
            echo -e "$accession\t$header\t$segment"
        done >> ${OUTPUT_DIR}/metadata/metadata.tsv
        
        echo "      âœ“ Simple metadata.tsv also created for compatibility"
        
        # Show summary
        total_records=$(tail -n +2 "${OUTPUT_DIR}/metadata/rich_metadata.csv" | wc -l)
        echo "      Rich metadata contains $total_records records with fields matching H5N1 CSV format"
        
        # Show sample isolate names
        echo "      Sample isolate names from metadata:"
        tail -n +2 "${OUTPUT_DIR}/metadata/rich_metadata.csv" | head -3 | cut -d',' -f1,10 | sed 's/^/        /'
    else
        echo "      âš ï¸  Rich metadata extraction failed, creating simple metadata only"
        
        # Fall back to simple metadata
        echo -e "accession\tfull_header\tsegment" > ${OUTPUT_DIR}/metadata/metadata.tsv
        
        grep "^>" "$SEQUENCES_FILE" | \
        while read header; do
            accession=$(echo "$header" | sed 's/^>//' | cut -d' ' -f1)
            segment=$(echo "$header" | grep -o "segment [0-9]" | grep -o "[0-9]" || echo "unknown")
            echo -e "$accession\t$header\t$segment"
        done >> ${OUTPUT_DIR}/metadata/metadata.tsv
    fi
else
    echo "      âš ï¸  JSON file not found, creating simple metadata only"
    
    # Create simple metadata as fallback
    echo -e "accession\tfull_header\tsegment" > ${OUTPUT_DIR}/metadata/metadata.tsv
    
    grep "^>" "$SEQUENCES_FILE" | \
    while read header; do
        accession=$(echo "$header" | sed 's/^>//' | cut -d' ' -f1)
        segment=$(echo "$header" | grep -o "segment [0-9]" | grep -o "[0-9]" || echo "unknown")
        echo -e "$accession\t$header\t$segment"
    done >> ${OUTPUT_DIR}/metadata/metadata.tsv
fi

# Cleanup
rm -f /tmp/extract_rich_metadata.py