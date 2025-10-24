# antimicrobial_collection

Code, pipelines and utilities for processing, featurizing, training and evaluating machine-learning models for antimicrobial/antibacterial screening datasets.

This repository collects several public screening datasets, standardizes SMILES, generates chemical representations (fingerprints), trains models, and evaluates results with reusable Nextflow components.

## General idea

- Inputs: CSV dataset(s) with at least a SMILES column and an activity label column (defaults controlled via `nextflow.config` — see `smiles_col` and `activity_col`).
- Outputs: per-run files in `results/` (for example `trained_model.joblib`, `results.csv`) and a collected summary `results/summary/all_results_summary.csv`.
- Errors: scripts exit non-zero on malformed input or missing dependencies; Nextflow reports process failures and logs.

## Repository layout

- `main.nf` — Orchestrating Nextflow workflow that composes modular components.
- `nextflow.config` — Default pipeline parameters (datasets path, component paths, columns, outdir, PYTHONPATH).
- `components/` — Component implementations grouped by stage:
  - `01_splitters/` — splitting strategies (random, scaffold, etc.)
  - `02_featurizers/` — featurizers (Morgan, MACCS, etc.)
  - `03_models/` — training scripts (e.g. `rf.py` for RandomForest)
  - `04_evaluators/` — evaluation metrics and reporting
- `src/` — shared Python helpers and data processing utilities (SMILES standardization, representations)
- `data/raw/` and `data/processed/` — datasets (processed CSVs are included for some public datasets; large raw datasets are omitted)
- `notebooks/` — example analysis and preprocessing notebooks
- `modules/` — Nextflow modular wrappers for running component groups (splitter/featurizer/trainer/evaluator)
- `results/` — example outputs from past runs (kept for reference)
- `DATA_SOURCES.md` — details and provenance for included datasets

## Requirements

To run the pipeline, you need nextflow and the conda environment.

```bash
git clone https://github.com/brunoholiva/antimicrobial_collection.git
cd antimicrobial_collection
conda env create -f env.yml
conda activate antimicrobial-collection
```

## Running the pipeline

From the repository root you can run the full pipeline with defaults in `nextflow.config`:

```bash
nextflow run main.nf
```

Nextflow will call the component scripts (splitters → featurizers → trainers → evaluators). Example per-process artifacts include `trained_model.joblib` and `results.csv`. At the end the workflow collects all evaluator results into `results/summary/all_results_summary.csv`.

Alternatively you can run changing specific parameters

```bash
nextflow run main.nf --outdir results --datasets ${PWD}/data/processed/dataset.csv --featurizers ${PWD}/components/02_featurizers/maccs.py
```

## Data sources

See `DATA_SOURCES.md` for dataset information and assay-level details for each included dataset. Large original datasets are not committed to this repo; processed CSVs for some public datasets are included under `data/processed/`.

## Outputs and results

- Per-run directories under `results/` or Nextflow work directories include trained models (`trained_model.joblib`) and per-run `results.csv` files.
- The workflow saves an aggregate summary to `results/summary/all_results_summary.csv` (see `main.nf` logic).

## Development and plans

This repository is in active development and will see changes.
Future plans include: including more models, featurizers, feature preprocessing, and plotting results.