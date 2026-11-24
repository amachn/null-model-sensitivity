# src/nullsens/backends/regioner.py
# - handles cross-language communication between Python and R, specifically regioneR

# libraries
import json
import subprocess
import sys

# project imports
from nullsens.analysis.models import AnalysisRequest, AnalysisResult
from nullsens.backends.runtime import find_rscript
from nullsens.paths import PACKAGE_PATHS


def run_regioner(request: AnalysisRequest) -> AnalysisResult:
    rscript = find_rscript()
    runner = PACKAGE_PATHS.r / "regioner_runner.R"

    payload = request.to_dict()

    # this method relies on us using "cat" for the actual output and "message" for diagnostics in R
    completed = subprocess.run(
        [rscript, runner], cwd=runner.parent, input=json.dumps(payload), capture_output=True, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"regioneR backend failed:\n{completed.stderr}")

    if request.diagnostics and completed.stderr:
        print(completed.stderr, file=sys.stderr, end="")

    result = json.loads(completed.stdout)
    return AnalysisResult(**result)
