# src/nullsens/analysis/runner.py
# - tool selection layer between the Python tool and the R backend

# project imports
from nullsens.analysis.models import AnalysisRequest, AnalysisResult
from nullsens.backends.regioner import run_regioner


def validate_request(request: AnalysisRequest):
    if not request.regions_a.is_file():
        raise FileNotFoundError(f"Region set A not found: {request.regions_a}")

    if not request.regions_b.is_file():
        raise FileNotFoundError(f"Region set B not found: {request.regions_b}")

    if not request.genome.strip():
        raise ValueError("Genome must be specified.")

    if request.permutations < 1:
        raise ValueError("Permutations must be greater than 0.")


def run_analysis(request: AnalysisRequest, backend: str = "regioner") -> AnalysisResult:
    validate_request(request)

    if backend == "regioner":
        return run_regioner(request)

    raise ValueError(f"Unsupported analysis backend: {backend}")
