# src/nullsens/analysis/models.py
# - defines the data classes that will manage requests and results

# libraries
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path


class LengthMode(StrEnum):
    EXACT = "exact"
    DISTRIBUTION = "distribution"


@dataclass
class NullConfig:
    length_mode: LengthMode = LengthMode.EXACT
    gc_matched: bool = False
    chromatin_matched: bool = False
    mappability_masked: bool = False


@dataclass
class ReferenceFiles:
    gc: Path | None = None
    chromatin: Path | None = None
    mappability: Path | None = None


@dataclass
class AnalysisRequest:
    regions_a: Path
    regions_b: Path
    genome: str

    permutations: int = 1000
    seed: int | None = None

    null: NullConfig = field(default_factory=NullConfig)
    references: ReferenceFiles = field(default_factory=ReferenceFiles)


@dataclass
class AnalysisResult:
    observed: float
    permutation_values: list[float]

    null_mean: float
    null_sd: float
    z_score: float
    p_value: float
