# antimicrobial_collection

This repository provides code and utilities for storing, processing, and analyzing molecules with antimicrobial properties, with a focus on high-throughput screening datasets.

## Features
- Preprocessing and standardization of SMILES strings
- Validation and deduplication of molecular datasets
- Utilities for handling multiple public and private antimicrobial screening datasets
- Example notebooks for data exploration

## Project Structure

- `src/` — Main source code (SMILES processing, utilities)
- `data/raw/` — Raw input datasets (not all included due to size)
- `data/processed/` — Processed datasets
- `misc/` — Helper scripts
- `preprocessing.ipynb` — Notebook with the actual preprocessing


## Data Sources

For detailed information about the datasets and their assay protocols, see [`DATA_SOURCES.md`](./DATA_SOURCES.md).