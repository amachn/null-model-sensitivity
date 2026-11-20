# src/nullsens/analysis/models.py
# - defines the data classes that will manage requests and results

# libraries
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ChromatinConfig:
    enabled: bool = False
    padding: int | None = None
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

