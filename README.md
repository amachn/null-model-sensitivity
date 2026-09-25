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
Large binary and genomic data files (like `.bed`, etc) are managed using Git LFS.

After initial setup, normal git commands can be used, and Git LFS will handle the tracked large files automatically.

#### Usage
To show the files currently managed w/ LFS:
```bash
git lfs ls-files
```

To show the file patterns currently tracked by LFS:
```bash
git lfs track
```

To add a new tracking pattern:
```bash
git lfs track "<pattern>"
```

For example:
```bash
git lfs track "*.bed"
git lfs track "*.bigwig"
git lfs track "*.xlsx"
```

A specific file/path can also be tracked:
```bash
git lfs track "reproductions/P03/results/*"
git lfs track "data/reference/encode_blacklist.bed"
```

Once a new pattern is tracked, you need to commit the `.gitattributes` file to update the repository.

To untrack a pattern, use the same syntax as above but with the `untrack` command: `git lfs untrack "<pattern>"`

## CLI
Current commands:

```bash
nullsens papers
 - add <paper-id [optional]>
 - edit <paper-id>
 - delete <paper-id>
 - show <paper-id>
 - list
 - refresh
```