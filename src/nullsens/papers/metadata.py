# src/nullsens/papers/metadata.py
# - handles the reading/writing of .yaml files used for paper metadata in root/papers

# libraries
from pathlib import Path
import re
import yaml

# project imports
from nullsens.paths import PATHS


# validate the provided paper ID
def normalize_paper_id(paper_id: str) -> str:
    paper_id = paper_id.upper().strip()

    if not re.fullmatch(r"P(?:0[1-9]|1[0-6])", paper_id):
        raise ValueError(f"Paper ID must be between P01 and P16.")

    return paper_id


# gets the papers directory
def get_paper_path(paper_id: str) -> Path:
    paper_id = normalize_paper_id(paper_id)
    return PATHS.papers / f"{paper_id.upper()}.yaml"


# fetch all existing paper YAMLs
def get_existing_paper_ids() -> list[str]:
    return [path.stem for path in sorted(PATHS.papers.glob("P*.yaml"))]


# gets the next paper ID in order from 01-16, accounting for gaps
# used to make the "nullsens papers add" paper_id arg optional
def get_next_paper_id() -> str:
    existing = set(get_existing_paper_ids())

    for number in range(1, 17):
        paper_id = f"P{number:02d}"

        if paper_id not in existing:
            return paper_id

    raise RuntimeError("All paper IDs P01-P16 are already in use.")


# loads the YAML file of the requested paper
def load_paper(paper_id: str) -> dict:
    path = get_paper_path(paper_id)

    if not path.exists():
        raise FileNotFoundError(f"Paper {paper_id} does not exist.")

    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# loads all papers, used to refresh the README file
def load_all_papers() -> list[dict]:
    papers = []

    for path in sorted(PATHS.papers.glob("P*.yaml")):
        with path.open("r", encoding="utf-8") as f:
            papers.append(yaml.safe_load(f))

    return papers


# updates the YAML file of the requested paper
def save_paper(paper_id: str, data: dict) -> str:
    paper_id = normalize_paper_id(paper_id)
    path = get_paper_path(paper_id)

    data["paper_id"] = paper_id
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)

    return paper_id


# creates a new data dictionary used in the creation of a paper's YAML file
def create_paper(paper_id: str) -> dict:
    paper_id = normalize_paper_id(paper_id)
    path = get_paper_path(paper_id)

    if path.exists():
        raise FileExistsError(f"Paper {paper_id} already exists.")

    return {
        "paper_id": paper_id,
        "title": "",
        "owner": "",
        "doi": "",
        "pmid": "",
        "data_url": "",
        "code_url": "",
        "reproduction_status": "Not Started",
    }


# deletes/drops a paper from the records
def drop_paper(paper_id: str) -> str:
    paper_id = normalize_paper_id(paper_id)
    path = get_paper_path(paper_id)

    if path.exists():
        path.unlink()
    else:
        raise FileNotFoundError(f"Paper {paper_id} does not exist.")

    return paper_id