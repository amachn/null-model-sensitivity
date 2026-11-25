# src/nullsens/analysis/reporting.py
# - formats, displays, and writes analysis results

# libraries
from pathlib import Path
import csv
import json
import math

# project imports
from nullsens.analysis.models import AnalysisRequest, AnalysisResult, FactorialCondition, FactorialResult

FACTORIAL_LABELS = {
    "baseline": "Baseline",
    "mappability": "Mappability",
    "gc": "GC",
    "chromatin": "Chromatin",
    "mappability_gc": "Mappability + GC",
    "mappability_chromatin": "Mappability + Chromatin",
    "gc_chromatin": "GC + Chromatin",
    "all_three": "All",
}

CSV_FIELDS = (
    "analysis_type",
    "null_model",
    "null_model_label",
    "status",
    "mappability",
    "gc",
    "chromatin",
    "n_regions",
    "observed",
    "expected",
    "null_sd",
    "log2_fold_change",
    "z_score",
    "p_value",
    "genome",
    "permutations",
    "seed",
    "cores",
    "max_attempts",
    "runtime_seconds",
    "error",
)

FACTORIAL_EFFECT_FIELDS = (
    "response",
    "term",
    "order",
    "effect",
    "absolute_effect",
)

FACTORIAL_EFFECT_TERMS = (
    ("mappability", ("mappability", )),
    ("gc", ("gc", )),
    ("chromatin", ("chromatin", )),
    ("mappability:gc", ("mappability", "gc")),
    ("mappability:chromatin", ("mappability", "chromatin")),
    ("gc:chromatin", ("gc", "chromatin")),
    ("mappability:gc:chromatin", ("mappability", "gc", "chromatin")),
)

# - output helpers -

# JSON doesn't support NaN/Inf when strict output is requested
def _clean_number(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None

    return value


# normalize the user-provided output prefix and create its parent directory
def _normalize_output_prefix(prefix: Path) -> Path:
    prefix = Path(prefix)

    if prefix.suffix.lower() in {".json", ".csv"}:
        prefix = prefix.with_suffix("")

    prefix.parent.mkdir(parents=True, exist_ok=True)

    return prefix


# derive the internal null-model key for a single analysis
def _null_model_key(mappability: bool, gc: bool, chromatin: bool) -> str:
    enabled = []

    if mappability:
        enabled.append("mappability")

    if gc:
        enabled.append("gc")

    if chromatin:
        enabled.append("chromatin")

    if not enabled:
        return "baseline"

    if len(enabled) == 3:
        return "all_three"

    return "_".join(enabled)


# convert one successful analysis result into JSON-safe data
def _analysis_result_dict(result: AnalysisResult) -> dict:
    return {
        "n_peaks": result.n_peaks,
        "observed": _clean_number(result.observed),
        "permutation_values": [_clean_number(value) for value in result.permutation_values],
        "null_mean": _clean_number(result.null_mean),
        "null_sd": _clean_number(result.null_sd),
        "log2_fold_change": _clean_number(result.log2_fold_change),
        "z_score": _clean_number(result.z_score),
        "p_value": _clean_number(result.p_value),
        "runtime_seconds": _clean_number(result.runtime_seconds),
    }


# construct one summary row for CSV output
def _csv_row(
    request: AnalysisRequest, *, analysis_type: str, null_model: str,
    mappability: bool, gc: bool, chromatin: bool,
    result: AnalysisResult | None, error: str | None = None
) -> dict:
    label = FACTORIAL_LABELS.get(null_model, null_model)

    return {
        "analysis_type": analysis_type,
        "null_model": null_model,
        "null_model_label": label,
        "status": "success" if result is not None else "failed",
        "mappability": mappability,
        "gc": gc,
        "chromatin": chromatin,
        "n_regions": result.n_peaks if result is not None else None,
        "observed": _clean_number(result.observed) if result is not None else None,
        "expected": _clean_number(result.null_mean) if result is not None else None,
        "null_sd": _clean_number(result.null_sd) if result is not None else None,
        "log2_fold_change": _clean_number(result.log2_fold_change) if result is not None else None,
        "z_score": _clean_number(result.z_score) if result is not None else None,
        "p_value": _clean_number(result.p_value) if result is not None else None,
        "genome": request.genome,
        "permutations": request.permutations,
        "seed": request.seed,
        "cores": request.cores,
        "max_attempts": request.max_attempts,
        "runtime_seconds": _clean_number(result.runtime_seconds) if result is not None else None,
        "error": error,
    }


def calculate_factorial_effects(factorial: FactorialResult) -> list[dict] | None:
    # all 8 successful cells are required for factorial contrasts
    if len(factorial.runs) != 8 or factorial.failed_runs:
        return None

    cells = []

    for run in factorial.runs:
        result = run.result

        if result is None:
            return None

        cells.append((run.condition, result.log2_fold_change))

    # ensure the full 2^3 design is represented exactly once
    conditions = {
        (condition.mappability, condition.gc, condition.chromatin) for condition, _ in cells
    }

    if len(conditions) != 8:
        return None

    effects = []

    for term, factors in FACTORIAL_EFFECT_TERMS:
        positive = []
        negative = []

        for condition, response in cells:
            levels = {
                "mappability": condition.mappability,
                "gc": condition.gc,
                "chromatin": condition.chromatin,
            }

            sign = 1
            for factor in factors:
                sign *= 1 if levels[factor] else -1

            if sign > 0:
                positive.append(response)
            else:
                negative.append(response)

        effect = (
                sum(positive) / len(positive)
                - sum(negative) / len(negative)
        )

        effects.append({
            "response": "log2_fold_change",
            "term": term,
            "order": len(factors),
            "effect": effect,
            "absolute_effect": abs(effect),
        })

    return effects


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")


def _write_factorial_effects_csv(path: Path, effects: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FACTORIAL_EFFECT_FIELDS)
        writer.writeheader()
        writer.writerows(effects)


# - public output writers -
def write_analysis_outputs(prefix: Path, request: AnalysisRequest, result: AnalysisResult) -> tuple[Path, Path]:
    prefix = _normalize_output_prefix(prefix)

    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")

    null_model = _null_model_key(
        request.null.mappability.enabled,
        request.null.gc.enabled,
        request.null.chromatin.enabled,
    )

    payload = {
        "analysis_type": "analyze",
        "request": request.to_dict(),
        "null_model": {
            "key": null_model,
            "label": FACTORIAL_LABELS.get(null_model, null_model),
        },
        "result": _analysis_result_dict(result),
    }

    row = _csv_row(
        request,
        analysis_type="analyze",
        null_model=null_model,
        mappability=request.null.mappability.enabled,
        gc=request.null.gc.enabled,
        chromatin=request.null.chromatin.enabled,
        result=result,
    )

    _write_json(json_path, payload)
    _write_csv(csv_path, [row])
    return json_path, csv_path


def write_factorial_outputs(
    prefix: Path,
    request: AnalysisRequest,
    factorial: FactorialResult
) -> tuple[Path, Path, Path | None]:
    prefix = _normalize_output_prefix(prefix)

    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    effects_path = prefix.with_suffix(".effects.csv")

    runs = []
    rows = []

    # build normal per-condition JSON + CSV data
    for run in factorial.runs:
        condition = run.condition

        runs.append({
            "condition": {
                "key": condition.label,
                "label": FACTORIAL_LABELS.get(condition.label, condition.label),
                "mappability": condition.mappability,
                "gc": condition.gc,
                "chromatin": condition.chromatin,
            },
            "status": "success" if run.succeeded else "failed",
            "error": run.error,
            "result": _analysis_result_dict(run.result) if run.result is not None else None,
        })

        rows.append(_csv_row(
            request,
            analysis_type="factorial",
            null_model=condition.label,
            mappability=condition.mappability,
            gc=condition.gc,
            chromatin=condition.chromatin,
            result=run.result,
            error=run.error,
        ))

    # calculate derived factorial contrasts
    effects = calculate_factorial_effects(factorial)

    # full archival JSON
    payload = {
        "analysis_type": "factorial",
        "request": request.to_dict(),
        "total_runtime_seconds": _clean_number(factorial.runtime_seconds),
        "effects": effects,
        "runs": runs,
    }

    # always write normal outputs
    _write_json(json_path, payload)
    _write_csv(csv_path, rows)

    # effects require a successful 2^3 factorial
    written_effects_path = None

    if effects is not None:
        _write_factorial_effects_csv(effects_path, effects)
        written_effects_path = effects_path

    return json_path, csv_path, written_effects_path


# - terminal formatting -
def format_runtime(result: AnalysisResult | FactorialResult) -> str:
    if result.runtime_seconds <= 60:
        return f"{result.runtime_seconds:.2f} s"

    return f"{result.runtime_seconds:.2f} s ({result.runtime_minutes:.2f} min)"


def print_analysis_result(result: AnalysisResult) -> None:
    print("Analysis complete.\n")
    print(f"Regions analyzed: {result.n_peaks}")
    print(f"Observed overlap: {result.observed:g}")
    print(f"Expected overlap: {result.null_mean:.4f}")
    print(f"Log2 fold change: {result.log2_fold_change:.4f}")
    print(f"Null SD: {result.null_sd:.4f}")
    print(f"Z-score: {result.z_score:.4f}")
    print(f"P-value: {result.p_value:.6g}")
    print(f"Runtime: {format_runtime(result)}")


def print_factorial_progress(idx: int, total: int, condition: FactorialCondition) -> None:
    label = FACTORIAL_LABELS.get(condition.label, condition.label)
    print(f"[{idx}/{total}] Running {label}...", flush=True)


def print_factorial_result(factorial: FactorialResult) -> None:
    print("Factorial analysis complete.\n")

    # table output header row
    headers = ("Null model", "N", "Observed", "Expected", "log2FC", "Z-score", "P-value", "Runtime")
    rows = []

    # iterate through each null result in the factorial
    for run in factorial.runs:
        label = FACTORIAL_LABELS.get(run.condition.label, run.condition.label)

        if not run.succeeded:
            rows.append((label, "-", "-", "-", "-", "-", "-", "FAILED"))
            continue

        result = run.result
        assert result is not None

        rows.append((
            label,
            str(result.n_peaks),
            f"{result.observed:g}",
            f"{result.null_mean:.2f}",
            f"{result.log2_fold_change:.4f}",
            f"{result.z_score:.2f}",
            f"{result.p_value:.6g}",
            format_runtime(result),
        ))

    # - table printing logic -

    # calculate width of each column, checking if the header or value is largest
    widths = [max(len(headers[idx]), *(len(row[idx]) for row in rows)) for idx in range(len(headers))]

    # output each row of the table left justified by the width of each column
    def format_row(row):
        return " ".join(value.ljust(widths[idx]) for idx, value in enumerate(row))

    print(format_row(headers))
    print("-" * sum(widths, len(headers) - 1))
    for row in rows:
        print(format_row(row))

    print(f"\nTotal analysis runtime: {format_runtime(factorial)}")

    if factorial.failed_runs:
        print("\nFailed runs:")

        for run in factorial.failed_runs:
            label = FACTORIAL_LABELS.get(run.condition.label, run.condition.label)
            error = run.error or "Unknown error"
            print(f"\n{label}: \n{error}")