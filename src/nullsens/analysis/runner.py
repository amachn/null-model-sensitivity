# src/nullsens/analysis/runner.py
# - tool selection layer between the Python tool and the R backend

# project imports
from nullsens.analysis.models import AnalysisRequest, AnalysisResult
from nullsens.backends.regioner import run_regioner
from nullsens.genomes.manager import is_genome_installed


def validate_request(request: AnalysisRequest):
    if not request.regions_a.is_file():
        raise FileNotFoundError(f"Region set A not found: {request.regions_a}")

    if not request.regions_b.is_file():
        raise FileNotFoundError(f"Region set B not found: {request.regions_b}")

    if not request.genome.strip():
        raise ValueError("Genome must be specified.")

    if not is_genome_installed(request.genome):
        raise RuntimeError(
            f"Genome '{request.genome}' is not installed.\n"
            f"Run: nullsens genomes install {request.genome}"
        )

    if request.permutations < 1:
        raise ValueError("Permutations must be greater than 0.")

    if request.max_attempts < 1:
        raise ValueError("Maximum attempts must be greater than 0.")

    if request.null.mappability.enabled:
        if request.references.mappability is None:
            raise ValueError("Mappability is enabled but no mappability reference was provided.")

        if not request.references.mappability.exists():
            raise ValueError(f"Mappability reference doesn't exist: {request.references.mappability}")

    if request.null.chromatin.enabled:
        if request.references.chromatin is None:
            raise ValueError("Chromatin matching is enabled but no chromatin reference was provided.")

        if not request.references.chromatin.exists():
            raise ValueError(f"Chromatin reference doesn't exist: {request.references.chromatin}")

    if request.null.gc.enabled:
        if request.null.gc.tolerance <= 0:
            raise ValueError("GC tolerance must be greater than 0.")

        if request.null.gc.flank < 0:
            raise ValueError("GC flank must be >= 0.")

        if not 0 <= request.null.gc.max_n <= 1:
            raise ValueError("GC max_n must be between 0 and 1.")


def run_analysis(request: AnalysisRequest, backend: str = "regioner") -> AnalysisResult:
    validate_request(request)

    if backend == "regioner":
        return run_regioner(request)

    raise ValueError(f"Unsupported analysis backend: {backend}")
