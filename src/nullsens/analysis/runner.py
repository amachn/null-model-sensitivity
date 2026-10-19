# src/nullsens/analysis/runner.py
# - tool selection layer between the Python tool and the R backend

# project imports
from nullsens.analysis.models import AnalysisRequest, AnalysisResult
from nullsens.backends import run_regioner


def run_analysis(request: AnalysisRequest, backend: str = "regioner") -> AnalysisResult:
    if backend == "regioner":
        return run_regioner(request)

    raise ValueError(f"Unsupported analysis backend: {backend}")
