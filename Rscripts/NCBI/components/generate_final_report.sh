#!/bin/bash

# Final Report Generation Component
# Part of WAPHL NCBI Pipeline

START_DATE="$1"

echo "Generating final comprehensive report..."

# Count sequences by subtype and segment
count_sequences() {
    local subtype=$1
    if [ -f "ncbi/$subtype/metadata/metadata.tsv" ]; then
        tail -n +2 "ncbi/$subtype/metadata/metadata.tsv" | wc -l
    else
        echo 0
    fi
}

# Get segment counts for a subtype
get_segment_counts() {
    local subtype=$1
    local segments=("pb2" "pb1" "pa" "ha" "np" "na" "mp" "ns")
    
    for segment in "${segments[@]}"; do
        if [ -f "ncbi/$subtype/sequences/${segment}_sequences.fasta" ]; then
            count=$(grep -c "^>" "ncbi/$subtype/sequences/${segment}_sequences.fasta" 2>/dev/null || echo 0)
            printf "  %-4s: %3d sequences\n" "$segment" "$count"
        else
            printf "  %-4s: %3d sequences\n" "$segment" "0"
        fi
    done
}

H3N2_TOTAL=$(count_sequences h3n2)
H1N1_TOTAL=$(count_sequences h1n1pdm)
VIC_TOTAL=$(count_sequences vic)
GRAND_TOTAL=$((H3N2_TOTAL + H1N1_TOTAL + VIC_TOTAL))

# Check QC status
QC_STATUS="Not run"
QC_ISSUES="Unknown"
if [ -f "qc_analysis/reports/qc_summary.txt" ]; then
    if grep -q "✅ PROCEED" qc_analysis/reports/qc_summary.txt; then
        QC_STATUS="✅ PASSED"
        QC_ISSUES="No critical issues"
    elif grep -q "❌ MANUAL REVIEW" qc_analysis/reports/qc_summary.txt; then
        QC_STATUS="❌ CRITICAL ISSUES"
        QC_ISSUES="Manual review required"
    elif grep -q "⚠️  REVIEW RECOMMENDED" qc_analysis/reports/qc_summary.txt; then
        QC_STATUS="⚠️  REVIEW NEEDED"
        QC_ISSUES="Issues detected - see QC report"
    else
        QC_STATUS="⚠️ WARNINGS"
        QC_ISSUES="Minor issues detected"
    fi
fi

# Generate comprehensive final report
cat > final_report.txt << EOF
================================================================
    WAPHL NCBI INFLUENZA PROCESSING PIPELINE - FINAL REPORT
================================================================
Date: $(date)
Start Date Filter: $START_DATE

PROCESSING SUMMARY:
------------------
Total WAPHL sequences processed: $GRAND_TOTAL

By Subtype:
- H3N2 sequences: $H3N2_TOTAL
- H1N1pdm sequences: $H1N1_TOTAL  
- Victoria sequences: $VIC_TOTAL

H3N2 SEGMENTS:
$(get_segment_counts h3n2)

H1N1PDM SEGMENTS:  
$(get_segment_counts h1n1pdm)

VICTORIA SEGMENTS:
$(get_segment_counts vic)

QUALITY CONTROL:
---------------
Status: $QC_STATUS
Issues: $QC_ISSUES
$(if [ -f "qc_analysis/reports/qc_summary.txt" ]; then
echo "Detailed QC Report: qc_analysis/reports/qc_summary.txt"
fi)

DIRECTORY STRUCTURE:
-------------------
ncbi/
├── h3n2/
│   ├── sequences/
│   │   ├── ha_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/ha_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── na_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/na_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── pb2_sequences.fasta     ($(grep -c "^>" ncbi/h3n2/sequences/pb2_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── pb1_sequences.fasta     ($(grep -c "^>" ncbi/h3n2/sequences/pb1_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── pa_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/pa_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── np_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/np_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   ├── mp_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/mp_sequences.fasta 2>/dev/null || echo 0) sequences)
│   │   └── ns_sequences.fasta      ($(grep -c "^>" ncbi/h3n2/sequences/ns_sequences.fasta 2>/dev/null || echo 0) sequences)
│   └── metadata/
│       ├── metadata.tsv            (sequence metadata)
│       └── nextclade_sort_results.tsv (classification results)
├── h1n1pdm/
│   ├── sequences/ (same structure as h3n2)
│   └── metadata/  (same structure as h3n2)
└── vic/
    ├── sequences/ (same structure as h3n2)
    └── metadata/
        └── metadata.tsv            (no NextClade - all Victoria)

qc_analysis/
├── reports/
│   └── qc_summary.txt              (comprehensive QC report)
└── flags/
    ├── cross_subtype_duplicates.txt (contamination issues)
    ├── low_confidence_sequences.txt (poor quality sequences)
    └── empty_segments.txt          (missing segments)

NEXT STEPS:
----------
$(if [ "$QC_STATUS" = "❌ CRITICAL ISSUES" ]; then
echo "🚨 CRITICAL ISSUES DETECTED:
1. ❌ DO NOT PROCEED with analysis until issues are resolved
2. Review QC report: qc_analysis/reports/qc_summary.txt
3. Investigate flagged sequences in qc_analysis/flags/
4. Address cross-contamination or processing errors
5. Re-run pipeline after corrections

CRITICAL ACTION REQUIRED BEFORE PROCEEDING"
elif [ "$QC_STATUS" = "⚠️  REVIEW NEEDED" ]; then
echo "⚠️  REVIEW RECOMMENDED:
1. Review QC report: qc_analysis/reports/qc_summary.txt
2. Consider excluding low-confidence sequences
3. After review, proceed with GISAID integration:
   - Merge NCBI sequences with GISAID by segment and subtype
   - Run Nextstrain builds for each subtype
   - Perform phylogenetic analysis"
else
echo "✅ READY TO PROCEED:
1. ✓ QC checks passed
2. Proceed with GISAID integration:
   - Merge NCBI sequences with GISAID by segment and subtype
   - Run Nextstrain builds for each subtype
   - Perform phylogenetic analysis"
fi)

KEY FILE LOCATIONS:
------------------
H3N2 Sequences:
  HA: ncbi/h3n2/sequences/ha_sequences.fasta
  NA: ncbi/h3n2/sequences/na_sequences.fasta
  Other segments: ncbi/h3n2/sequences/[segment]_sequences.fasta
  Metadata: ncbi/h3n2/metadata/metadata.tsv

H1N1pdm Sequences:
  HA: ncbi/h1n1pdm/sequences/ha_sequences.fasta  
  NA: ncbi/h1n1pdm/sequences/na_sequences.fasta
  Other segments: ncbi/h1n1pdm/sequences/[segment]_sequences.fasta
  Metadata: ncbi/h1n1pdm/metadata/metadata.tsv

Victoria Sequences:
  HA: ncbi/vic/sequences/ha_sequences.fasta
  NA: ncbi/vic/sequences/na_sequences.fasta
  Other segments: ncbi/vic/sequences/[segment]_sequences.fasta
  Metadata: ncbi/vic/metadata/metadata.tsv

PIPELINE STATUS: $(if [ "$QC_STATUS" = "✅ PASSED" ]; then echo "✅ COMPLETED SUCCESSFULLY"; elif [ "$QC_STATUS" = "❌ CRITICAL ISSUES" ]; then echo "❌ COMPLETED WITH CRITICAL ISSUES"; else echo "⚠️ COMPLETED WITH WARNINGS"; fi)

================================================================
Pipeline completed: $(date)
================================================================
EOF

echo "✓ Final report generated: final_report.txt"

# Also create a simple summary for quick reference
cat > processing_summary.txt << EOF
WAPHL NCBI Processing Summary
============================
Date: $(date)
Total sequences: $GRAND_TOTAL
- H3N2: $H3N2_TOTAL
- H1N1pdm: $H1N1_TOTAL  
- Victoria: $VIC_TOTAL
QC Status: $QC_STATUS
$(if [ "$QC_STATUS" != "✅ PASSED" ]; then
echo "⚠️  Review required: qc_analysis/reports/qc_summary.txt"
fi)
EOF

echo "✓ Summary generated: processing_summary.txt"
