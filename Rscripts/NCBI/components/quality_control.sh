#!/bin/bash

# Quality Control Component
# Part of WAPHL NCBI Pipeline

set -e

echo "Running quality control analysis..."

# Create QC directory
mkdir -p qc_analysis/{reports,flags}

# Check for cross-subtype contamination
echo "  Checking for cross-subtype contamination..."
DUPLICATES=0
if [ -f "ncbi/h3n2/metadata/metadata.tsv" ] && [ -f "ncbi/h1n1pdm/metadata/metadata.tsv" ]; then
    # Extract accessions from each subtype
    tail -n +2 ncbi/h3n2/metadata/metadata.tsv | cut -f1 | sort > qc_analysis/h3n2_accessions.txt
    tail -n +2 ncbi/h1n1pdm/metadata/metadata.tsv | cut -f1 | sort > qc_analysis/h1n1_accessions.txt
    
    # Find duplicates and create informative output
    comm -12 qc_analysis/h3n2_accessions.txt qc_analysis/h1n1_accessions.txt > qc_analysis/flags/cross_subtype_duplicates.txt
    
    DUPLICATES=$(wc -l < qc_analysis/flags/cross_subtype_duplicates.txt)
    if [ "$DUPLICATES" -gt 0 ]; then
        echo "    ⚠️  WARNING: Found $DUPLICATES sequences in both H3N2 and H1N1pdm data!"
        echo "    See: qc_analysis/flags/cross_subtype_duplicates.txt"
        echo "    Duplicated accessions:"
        head -5 qc_analysis/flags/cross_subtype_duplicates.txt | sed 's/^/      /'
        if [ "$DUPLICATES" -gt 5 ]; then
            echo "      ... and $((DUPLICATES - 5)) more"
        fi
    else
        echo "    ✓ No cross-subtype contamination detected"
    fi
else
    echo "    ⚠️  Missing metadata files - skipping cross-contamination check"
fi

# Check NextClade sort confidence if available
echo "  Analyzing NextClade confidence scores..."
LOW_CONFIDENCE=0
if [ -f "temp_processing/sort_results.tsv" ]; then
    # Create header for low confidence file
    echo -e "accession\tdataset\tconfidence_score" > qc_analysis/flags/low_confidence_sequences.txt
    
    # Extract low confidence sequences (score < 0.8) with accession numbers
    awk -F'\t' 'NR>1 && $4 < 0.8 {
        # Extract accession from seqName (column 2) - everything before the first space
        split($2, parts, " ");
        accession = parts[1];
        print accession "\t" $3 "\t" $4
    }' temp_processing/sort_results.tsv >> qc_analysis/flags/low_confidence_sequences.txt
    
    LOW_CONFIDENCE=$(tail -n +2 qc_analysis/flags/low_confidence_sequences.txt | wc -l)
    
    if [ "$LOW_CONFIDENCE" -gt 0 ]; then
        echo "    ⚠️  Found $LOW_CONFIDENCE sequences with low confidence scores"
        echo "    See: qc_analysis/flags/low_confidence_sequences.txt"
        echo "    Sample low confidence sequences:"
        head -5 qc_analysis/flags/low_confidence_sequences.txt | tail -n +2 | sed 's/^/      /'
    else
        echo "    ✓ All sequences have good confidence scores"
    fi
else
    echo "    ⚠️  NextClade sort results not found - skipping confidence analysis"
fi

# Check for empty segment files
echo "  Checking segment file completeness..."
EMPTY_SEGMENTS=0
> qc_analysis/flags/empty_segments.txt
for subtype in h3n2 h1n1pdm vic; do
    if [ -d "ncbi/$subtype/sequences" ]; then
        for segment_file in ncbi/$subtype/sequences/*.fasta; do
            if [ -f "$segment_file" ]; then
                count=$(grep -c "^>" "$segment_file" 2>/dev/null || echo 0)
                if [ "$count" -eq 0 ]; then
                    echo "$subtype/$(basename $segment_file .fasta)" >> qc_analysis/flags/empty_segments.txt
                    EMPTY_SEGMENTS=$((EMPTY_SEGMENTS + 1))
                fi
            fi
        done
    fi
done

if [ "$EMPTY_SEGMENTS" -gt 0 ]; then
    echo "    ⚠️  Found $EMPTY_SEGMENTS empty segment files"
    echo "    See: qc_analysis/flags/empty_segments.txt"
else
    echo "    ✓ All segment files contain sequences"
fi

# Count total sequences processed
H3N2_SEQS=0
H1N1_SEQS=0
VIC_SEQS=0

if [ -f "ncbi/h3n2/metadata/metadata.tsv" ]; then
    H3N2_SEQS=$(tail -n +2 ncbi/h3n2/metadata/metadata.tsv | wc -l)
fi

if [ -f "ncbi/h1n1pdm/metadata/metadata.tsv" ]; then
    H1N1_SEQS=$(tail -n +2 ncbi/h1n1pdm/metadata/metadata.tsv | wc -l)
fi

if [ -f "ncbi/vic/metadata/metadata.tsv" ]; then
    VIC_SEQS=$(tail -n +2 ncbi/vic/metadata/metadata.tsv | wc -l)
fi

# Generate QC summary report
cat > qc_analysis/reports/qc_summary.txt << EOF
WAPHL NCBI Pipeline - Quality Control Report
===========================================
Date: $(date)

SEQUENCES PROCESSED:
- H3N2 sequences: $H3N2_SEQS
- H1N1pdm sequences: $H1N1_SEQS  
- Victoria sequences: $VIC_SEQS
- Total sequences: $((H3N2_SEQS + H1N1_SEQS + VIC_SEQS))

QUALITY CONTROL RESULTS:
- Cross-subtype duplicates: $DUPLICATES sequences
- Low confidence sequences: $LOW_CONFIDENCE sequences  
- Empty segment files: $EMPTY_SEGMENTS files

$(if [ "$DUPLICATES" -gt 0 ]; then
echo "🚨 CRITICAL ISSUES FOUND:
- Sequences appear in multiple subtypes (see qc_analysis/flags/cross_subtype_duplicates.txt)
- Manual review required before proceeding with analysis
- These sequences may indicate:
  * Cross-contamination during sequencing
  * Reassortant viruses with mixed segments
  * Processing errors in classification"
fi)

$(if [ "$LOW_CONFIDENCE" -gt 0 ]; then
echo "⚠️  LOW CONFIDENCE SEQUENCES:
- Some sequences have NextClade confidence scores < 0.8
- Review qc_analysis/flags/low_confidence_sequences.txt
- Consider excluding from phylogenetic analysis"
fi)

$(if [ "$EMPTY_SEGMENTS" -gt 0 ]; then
echo "📋 MISSING SEGMENTS:
- Some expected segments have no sequences
- See qc_analysis/flags/empty_segments.txt
- May indicate limited sequencing coverage"
fi)

RECOMMENDATION:
$(if [ "$DUPLICATES" -gt 0 ]; then
echo "❌ MANUAL REVIEW REQUIRED - Critical cross-contamination detected"
elif [ "$LOW_CONFIDENCE" -gt 10 ]; then
echo "⚠️  REVIEW RECOMMENDED - Many low confidence sequences"
else
echo "✅ PROCEED - No critical issues detected"
fi)

QC FILES GENERATED:
- This report: qc_analysis/reports/qc_summary.txt
- Cross-subtype duplicates: qc_analysis/flags/cross_subtype_duplicates.txt
- Low confidence sequences: qc_analysis/flags/low_confidence_sequences.txt
- Empty segments: qc_analysis/flags/empty_segments.txt

$(if [ "$DUPLICATES" -gt 0 ] || [ "$LOW_CONFIDENCE" -gt 10 ]; then
echo "
NEXT STEPS:
1. Review all flagged sequences
2. Investigate causes of issues
3. Re-process or exclude problematic sequences
4. Re-run QC after corrections"
else
echo "
NEXT STEPS:
1. Proceed with GISAID integration
2. Run Nextstrain builds
3. Perform phylogenetic analysis"
fi)
EOF

echo "✓ Quality control analysis complete"

# Return appropriate exit code
if [ "$DUPLICATES" -gt 0 ]; then
    echo "  ⚠️  Critical QC issues found - manual review recommended"
    exit 1
else
    echo "  ✓ QC passed - no critical issues detected"
    exit 0
fi
