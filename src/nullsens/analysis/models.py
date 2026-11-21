# src/nullsens/analysis/models.py
# - defines the data classes that will manage requests and results

# libraries
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ChromatinConfig:
    enabled: bool = False
    padding: int = 0
    label_col: str | None = None
    priority: list[str] | None = None
    strict: bool = False


@dataclass
class GCConfig:
    enabled: bool = False
    tolerance: float = 0.02
    flank: int = 0
    max_n: float = 0.1


@dataclass
class MappabilityConfig:
    enabled: bool = False
    threshold: float | None = None

    # False -> reference contains forbidden/excluded regions
    # True -> reference contains allowed workspace; its complement becomes the forbidden mask
    inclusion_workspace: bool = False


@dataclass
class NullConfig:
    chromatin: ChromatinConfig = field(default_factory=ChromatinConfig)
    gc: GCConfig = field(default_factory=GCConfig)
    mappability: MappabilityConfig = field(default_factory=MappabilityConfig)


@dataclass
class ReferenceFiles:
    chromatin: Path | None = None
    mappability: Path | None = None


@dataclass
class AnalysisRequest:
    regions_a: Path
    regions_b: Path
    genome: str

    cores: int = 1
    permutations: int = 1000
    max_attempts: int = 100
    seed: int | None = None

    null: NullConfig = field(default_factory=NullConfig)
    references: ReferenceFiles = field(default_factory=ReferenceFiles)


@dataclass
class AnalysisResult:
    n_peaks: int
    observed: float
    permutation_values: list[float]

    null_mean: float
    null_sd: float

    log2_fold_change: float
    z_score: float
    p_value: float

    runtime_seconds: float

    @property
    def runtime_minutes(self) -> float:
        return self.runtime_seconds / 60.0


@dataclass(frozen=True)
class FactorialCondition:
    label: str
    mappability: bool
    gc: bool
    chromatin: bool


@dataclass
class FactorialRunResult:
    condition: FactorialCondition
    result: AnalysisResult | None = None
    error: str | None = None

    @property
    def succeeded(self) -> bool:
        return self.result is not None


@dataclass
class FactorialResult:
    runs: list[FactorialRunResult]

    @property
    def successful_runs(self) -> list[FactorialRunResult]:
        return [run for run in self.runs if run.succeeded]

    @property
    def failed_runs(self) -> list[FactorialRunResult]:
        return [run for run in self.runs if not run.succeeded]

    @property
    def runtime_seconds(self) -> float:
        return sum(
            run.result.runtime_seconds
            for run in self.successful_runs
            if run.result is not None
        )

    @property
    def runtime_minutes(self) -> float:
        return self.runtime_seconds / 60.0
