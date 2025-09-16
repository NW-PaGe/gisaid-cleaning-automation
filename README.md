# WA DOH GISAID Database Download and Cleaning for Flu Nextstrain Builds
This repository contains scripts for cleaning, merging, and maintaining the WADOH MEP flu database for avian flu and seasonal flu. This database pulls in data from NCBI and GISAID.
This repository is to be used in conjunction with our AWS resources.

## Features of the database
- Merges NCBI and GISAID Database
- Cleans HA segment for h5n1
- Cleans HA and NA segments for h3n2, h1n1pdm, Victoria


## Requirements:
- Linux based system
- conda installed, if not consider installing [miniconda](https://docs.anaconda.com/miniconda/install/)
- awscli

## GISAID Search Parameters for Each Influenza Type:

|               | h5n1          | h3n2 | h1n1pdm | victoria  |
|--------------|--------------|------|---------|-----------|
| **Type**     | A            | A    | A       | B         |
| **H**        | 5            | 3    | 1       | -         |
| **N**        | 1            | 2    | 1       | -         |
| **Lineage**  | -            | -    | pdm09   | Victoria  |
| **Host**     | Human, Animal | all  | all     | all       |
| **Location** | all          | all  | all     | all       |
| **Submission date from** | last 2 months | last 2 months | last 2 months | last 2 months |
| **Required segments** | HA  | -    | -       | -         |


### GISAID Sequence Download Steps
1. With the search parameters, download Sequences as FASTA with format
- Format: "Sequences (DNA) as FASTA"
- FASTA Header: Isolate ID
- Date format: YYYY-MM-DD
- Check replace spaces with underscores in FASTA header
- Check remove space before and after values in FASTA header
- For H5N1 Download HA segment only. For H3N2, H1N1, and Vic Download HA and NA segments into separate files.

2. After downloading sequences, save FASTA file in `data/raw/[flu-strain]/fasta/` folder as `[Year-present]_[Month-present]_[Day-present]_[Year-past]_[Month-past]_[Day-past]_[flu-strain]_[segment]_sequences`


### GISAID Metadata Download
1. With the same search parameters as the previous step, download "Isolates as XLS (virus metadata only)"
- Date format : YYYY-MM-DD

2. After downloading, save XLS file in `data/raw/[flu-strain]/metadata/` folder as `[Year-present]_[Month-present]_[Day-present]_[Year-past]_[Month-past]_[Day-past]_[flu-strain]_metadata`



## How to run the cleaning scripts
Once you have saved your GISAID metadata and fasta files into the respective flu-strain folders under `data/` and you've cloned this repo you'll move into this repository and run the `run_cleaning_gisaid.sh` script.
```
git clone <repo>

cd gisaid-database

chmod +x 01_run_cleaning_gisaid.sh

aws s3 sync [s3 bucket for data/raw/] data/raw/

./01_run_cleaning_gisaid.sh
```

### NCBI Data Downloads
- for sequence files, download through NCBI Virus
  - Step 1: Select Data Type
    - Sequence Data (FASTA format): nucleotide
  - Step 2: Select Records
  - Step 3: Select FASTA definition line
    - Use default: Accession GenBank Title
- For metadata files:
  - Step 1: Select Data Type
    - Results Table: CSV format
  - Step 2: Select Records
    - Download Selected Records or Download All Records
  - Step 3: Select columns to include in results Set
    - Select All
    - Accession: with version

## File Structure Overview
Pauline TODO: add in file tree


## Expected Output
Cleaned data will be outputted into `data/cleaned/[flu-strain]/fasta` and `data/cleaned/[flu-strain]/metadata`. The fasta folder will contain fasta files for ha and na sequences separately and the metadata file will contain the metadata matching all those sequences.

These data files are now ready to be ingested into the seasonal-flu and avian-flu nextstrain builds.

One additional folder can be found under `data/cleaned/gisaid_catalog` that will contain the combined catalog for all strains processed.  
