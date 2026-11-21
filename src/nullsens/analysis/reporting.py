# src/nullsens/analysis/reporting.py
# - formats and displays analysis results

# project imports
from nullsens.analysis.models import AnalysisResult, FactorialCondition, FactorialResult

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


def format_runtime(result: AnalysisResult | FactorialResult) -> str:
    if result.runtime_seconds <= 60:
        return f"{result.runtime_seconds:.2f} s"

    return f"{result.runtime_seconds:.2f} s ({result.runtime_minutes:.2f}"


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
            print(f"\n{label}: \n{run.error}")