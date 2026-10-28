# src/nullsens/backends/regioner.py
# - handles cross-language communication between Python and R, specifically regioneR

# libraries
from dataclasses import asdict
from pathlib import Path
import json
import shutil
import subprocess

# project imports
from nullsens.analysis import AnalysisRequest, AnalysisResult
from nullsens.paths import PACKAGE_PATHS


def run_regioner(request: AnalysisRequest) -> AnalysisResult:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found on PATH.")

    runner = PACKAGE_PATHS.r / "regioner_runner.R"

    payload = {
        "regions_a": str(request.regions_a.resolve()),
        "regions_b": str(request.regions_b.resolve()),
        "genome": request.genome,
        "permutations": request.permutations,
        "seed": request.seed,
        "null": {
            "length_mode": request.null.length_mode.value,
            "gc_matched": request.null.gc_matched,
            "chromatin_matched": request.null.chromatin_matched,
            "mappability_masked": request.null.mappability_masked,
        },
        "references": {
            name: (str(path.resolve()) if path else None)
            for name, path in asdict(request.references).items()
        }
    }

    # this method relies on us using "cat" for the actual output and "message" for diagnostics in R
    completed = subprocess.run(
        [rscript, runner], cwd=runner.parent, input=json.dumps(payload), capture_output=True, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"regioneR backend failed:\n{completed.stderr}")

    result = json.loads(completed.stdout)
    return AnalysisResult(**result)
