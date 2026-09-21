# Null Model Sensitivity

Assessment of null model sensitivity in published genomic enrichment analyses.

## Repository Structure

- `src/` = Python source code and CLI
- `papers/` = metadata for selected papers
- `reproductions/` = reproduction work for individual papers
- `data/` = shared/reference datasets (genomes, etc)

## Development Setup

### Requirements
- Git
- Git LFS
- Python >= 3.11

### Clone the Repository
Install and initalize Git LFS before cloning:
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

### Git LFS
Large data files (like `.bed`, etc) are managed using Git LFS.

After initial setup, normal git commands can be used, and Git LFS will handle tracked large files automatically.

## CLI
Current commands:

```bash
nullsens papers
 - add <paper-id [optional]>
 - edit <paper-id>
 - delete <paper-id>
 - list
 - refresh
```