# src/nullsens/backends/r_runtime.py
# - shared utilities for interacting with the local R install

# libraries
from pathlib import Path
import shutil


# simple shutil 'which' utility for finding the Rscript executable
def find_rscript() -> Path:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found on PATH.")
    return Path(rscript)
