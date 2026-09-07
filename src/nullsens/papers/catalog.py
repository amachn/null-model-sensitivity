# src/nullsens/papers/catalog.py
# - handles the status updating of the README.md file in root/papers

# libraries
from enum import Enum

# project imports
from nullsens.paths import PATHS
from nullsens.papers.metadata import load_all_papers

# statics
README = PATHS.papers / "README.md"
START_MARKER = "<!-- PAPER_TABLE_START -->"
END_MARKER = "<!-- PAPER_TABLE_END -->"


# used for hyperlink generation
class _LinkType(Enum):
    DOI = 1
    PMID = 2
    GENERIC = 3


# creates a specific hyperlink for github-flavored markdown
def _create_hyperlink(value: str, link: _LinkType) -> str:
    if not value:
        return ""

    value = value.strip()
    match link:
        case _LinkType.DOI:
            url = f"https://doi.org/{value}"
        case _LinkType.PMID:
            url = f"https://pubmed.ncbi.nlm.nih.gov/{value}"
        case _LinkType.GENERIC:
            url = value

    return f"[Link]({url})"


# builds the papers/README.md table
def build_table(papers: list[dict]) -> str:
    lines = [
        "| ID | Paper | Owner | DOI | PMID | Data | Code | Reproduction |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for paper in papers:
        paper_id = paper.get("paper_id", "")
        title = paper.get("title", "")
        owner = paper.get("owner", "")

        doi = _create_hyperlink(paper.get("doi", ""), _LinkType.DOI)
        pmid = _create_hyperlink(paper.get("pmid", ""), _LinkType.PMID)
        data = _create_hyperlink(paper.get("data_url", ""), _LinkType.GENERIC)
        code = _create_hyperlink(paper.get("code_url", ""), _LinkType.GENERIC)

        repro_status = paper.get("reproduction_status", "")

        lines.append(
            f"| {paper_id} | {title} | {owner} | {doi} |"
            f"{pmid} | {data} | {code} | {repro_status} |"
        )

    return "\n".join(lines)


# performs the actual README file refresh operation
def refresh_readme() -> None:
    papers = load_all_papers()
    table = build_table(papers)

    text = README.read_text(encoding="utf-8")

    if START_MARKER not in text or END_MARKER not in text:
        raise RuntimeError("Table markers are missing from papers/README.md.")

    before, remainder = text.split(START_MARKER, 1)
    _, after = remainder.split(END_MARKER, 1)

    updated = (
        before + START_MARKER + "\n\n"
        + table
        + "\n\n" + END_MARKER + after
    )

    README.write_text(updated, encoding="utf-8")
