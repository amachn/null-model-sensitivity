# src/nullsens/paths.py
# - this script makes resolving the project root + path aggregation simple

# libraries
from dataclasses import dataclass
from pathlib import Path

# blank file used to demarcate the top-level of this project
PROJECT_HEAD = ".capstone"


# performs an upward directory search from the run location to find the top directory
def find_project_root() -> Path:
    path = Path.cwd().resolve()
    if path.is_file():
        path = path.parent

    for candidate in (path, *path.parents):
        if (candidate / PROJECT_HEAD).exists():
            return candidate

    raise RuntimeError(f"Could not find project root ({PROJECT_HEAD}).")


# dataclass used to easily access the pathlib Paths for specific directories
@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @property
    def papers(self) -> Path:
        return self.root / "papers"

    @property
    def reproductions(self) -> Path:
        return self.root / "reproductions"

    @property
    def data(self) -> Path:
        return self.root / "data"

    # add future commonly used project paths (files too) here in the same manner as above


# pre-initialized ProjectPaths
PATHS = ProjectPaths(find_project_root())