#!/bin/bash

# WAPHL NCBI Influenza Processing Pipeline - Main Orchestrator
# Author: Pauline (WA DOH MEP)
# Purpose: Complete pipeline for downloading, processing, and QCing WAPHL influenza sequences

set -e

echo "============================================================"
echo "           WAPHL NCBI Influenza Processing Pipeline        "
echo "============================================================"
echo "Date: $(date)"
echo ""

# Configuration
START_DATE="2025-01-01"
OUTPUT_BASE="ncbi"
SKIP_QC=false
SKIP_INFLUENZA_B=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --start-date)
            START_DATE="$2"
            shift 2
            ;;
        --skip-qc)
            SKIP_QC=true
            shift
            ;;
        --skip-flu-b)
            SKIP_INFLUENZA_B=true
            shift
            ;;
        --help)
            cat << EOF
WAPHL NCBI Influenza Processing Pipeline

Usage: $0 [OPTIONS]

Options:
    --start-date YYYY-MM-DD    Start date for sequence filtering (default: 2025-01-01)
    --skip-qc                  Skip quality control analysis
    --skip-flu-b               Skip Influenza B processing
    --help                     Show this help message

Description:
    Complete pipeline for processing WAPHL influenza sequences from NCBI:
    1. Downloads and processes Influenza A (H3N2, H1N1pdm) 
    2. Downloads and processes Influenza B (Victoria)
    3. Runs quality control analysis
    4. Generates comprehensive reports

Output:
    ncbi/h3n2/sequences/        - H3N2 sequences by segment
    ncbi/h1n1pdm/sequences/     - H1N1pdm sequences by segment  
    ncbi/vic/sequences/         - Victoria sequences by segment
    qc_analysis/                - Quality control reports

EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "Configuration:"
echo "  Start date: $START_DATE"
echo "  Skip QC: $SKIP_QC"
echo "  Skip Influenza B: $SKIP_INFLUENZA_B"
echo ""

# Check dependencies
echo "Checking dependencies..."
required_tools=("datasets" "seqtk" "nextclade")
missing_tools=()

for tool in "${required_tools[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        missing_tools+=("$tool")
    fi
done

if [ ${#missing_tools[@]} -ne 0 ]; then
    echo "❌ Missing required tools: ${missing_tools[*]}"
    echo "Please install missing tools and try again."
    exit 1
fi

echo "✓ All required tools found"
echo ""

# Clean up any existing results
echo "Cleaning up previous results..."
rm -rf ${OUTPUT_BASE}
rm -rf qc_analysis
rm -rf temp_processing

echo ""
echo "=========================================="
echo "Phase 1: Processing Influenza A sequences"
echo "=========================================="

# Process Influenza A
echo "Starting Influenza A (H3N2/H1N1pdm) processing..."
if ! bash components/process_influenza_a.sh "$START_DATE"; then
    echo "❌ Influenza A processing failed"
    exit 1
fi

echo "✓ Influenza A processing completed"

# Process Influenza B
if [ "$SKIP_INFLUENZA_B" = false ]; then
    echo ""
    echo "=========================================="
    echo "Phase 2: Processing Influenza B sequences"
    echo "=========================================="
    
    echo "Starting Influenza B (Victoria) processing..."
    if ! bash components/process_influenza_b.sh "$START_DATE"; then
        echo "❌ Influenza B processing failed"
        exit 1
    fi
    
    echo "✓ Influenza B processing completed"
else
    echo ""
    echo "Skipping Influenza B processing (--skip-flu-b)"
fi

# Quality Control
if [ "$SKIP_QC" = false ]; then
    echo ""
    echo "========================================"
    echo "Phase 3: Quality Control Analysis"
    echo "========================================"
    
    echo "Running quality control analysis..."
    if ! bash components/quality_control.sh; then
        echo "⚠️  Quality control analysis failed (continuing anyway)"
    else
        echo "✓ Quality control analysis completed"
    fi
else
    echo ""
    echo "Skipping quality control analysis (--skip-qc)"
fi

echo ""
echo "========================================"
echo "Phase 4: Final Report Generation"
echo "========================================"

# Generate comprehensive final report
bash components/generate_final_report.sh "$START_DATE"

echo ""
echo "============================================================"
echo "                   PIPELINE COMPLETE!                      "
echo "============================================================"
echo ""
echo "Results organized in:"
echo "  ${OUTPUT_BASE}/              - All processed sequences"
echo "  qc_analysis/                 - Quality control reports"
echo "  final_report.txt             - Comprehensive summary"
echo ""
echo "Next steps:"
echo "  1. Review QC reports for any flagged sequences"
echo "  2. Merge with GISAID sequences by segment and subtype"
echo "  3. Run Nextstrain builds"
echo ""
echo "Pipeline completed successfully! 🎉"
