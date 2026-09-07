# Null Model Sensitivity

Assessment of null model sensitivity in published genomic enrichment analyses.

## Repository Structure

- `src/` = Python source code and CLI
- `papers/` = metadata for selected papers
- `reproductions/` = reproduction work for individual papers
- `data/` = shared/reference datasets (genomes, etc)

## Development Setup

Clone the repository, then create a virtual environment:

```bash
python -m venv .venv # python may be replaced by py or python3 depending on your setup
```

Activate the environment.

### Windows

```bash
.venv/Scripts/activate
```

### macOS/Linux

```bash
source .venv/bin/activate
```

Then, install the project in editable mode:

```bash
python -m pip install -e .
```

This installs the `nullsens` CLI while keeping the source code editable.

## CLI

Current commands:

```bash
nullsens papers
 - list
 - add <paper-id [optional]>
 - edit <paper-id>
 - refresh
```
