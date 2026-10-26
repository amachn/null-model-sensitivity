# Null Model Sensitivity

Assessment of null model sensitivity in published genomic enrichment analyses.

## Repository Structure

- `src/nullsens/` = Python source code and CLI
  - `analysis/` = shared analysis models and runners
  - `backends/` = Python interfaces to analysis backends
  - `papers/` = paper metadata management
  - `r/` = R backend, null-model implementations, and `renv` environment
- `papers/` = metadata for selected papers
- `reproductions/` = reproduction work for individual papers
- `data/` = shared/reference datasets (genomes, etc)
- `results/` = shared analysis and sensitivity results

## Installation

### Requirements
- Git
- Git LFS
- Python >= 3.11
- R
  - The current environment was created with R 4.6.1 and Bioconductor 3.23.
- `Rscript` available on your system PATH

### Clone the Repository
Install and initialize Git LFS before cloning:
```bash
git lfs install
```

Clone the repository:
```bash
git clone https://github.com/amachn/null-model-sensitivity.git
cd null-model-sensitivity
```

Ensure all Git LFS files are downloaded:
```bash
git lfs pull
```

### Python Environment
Create a virtual environment:
```bash
python -m venv .venv # python may be replaced by py or python3 depending on your setup
```

Then, activate the environment:
#### Windows
```bash
.venv\Scripts\activate
```

#### macOS/Linux
```bash
source .venv/bin/activate
```

Then, install the project in editable mode:

```bash
python -m pip install -e .
```

This installs the `nullsens` CLI while keeping the source code editable.

### R Environment
R dependencies are managed using `renv`. 

From the repository root, run:
```
cd src/nullsens/r
Rscript -e "renv::restore(prompt = FALSE)
```

### Git LFS
Large binary and genomic data files (like `.bed`, etc) are managed using Git LFS.

After initial setup, normal git commands can be used, and Git LFS will handle the tracked large files automatically.

## CLI
Use:

```
nullsens --help
```

to view available commands.