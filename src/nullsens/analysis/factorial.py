# src/nullsens/analysis/factorial.py
# - generates and runs factorial null-model configurations

# libraries
from dataclasses import replace

# project imports
from nullsens.analysis.models import (
    AnalysisRequest,
    FactorialCondition,
    FactorialResult,
    FactorialRunResult,
    NullConfig,
)
from nullsens.analysis.runner import run_analysis

# full 2^3 factorial design: baseline -> individual -> pairs -> all
FACTORIAL_DESIGN = (
    FactorialCondition("baseline",              mappability=False, gc=False, chromatin=False),
    FactorialCondition("mappability",           mappability=True,  gc=False, chromatin=False),
    FactorialCondition("gc",                    mappability=False, gc=True,  chromatin=False),
    FactorialCondition("chromatin",             mappability=False, gc=False, chromatin=True),
    FactorialCondition("mappability_gc",        mappability=True,  gc=True,  chromatin=False),
    FactorialCondition("mappability_chromatin", mappability=True,  gc=False, chromatin=True),
    FactorialCondition("gc_chromatin",          mappability=False, gc=True,  chromatin=True),
    FactorialCondition("all_three",             mappability=True,  gc=True,  chromatin=True),
)


# build a NullConfig matching the current factorial condition
def build_factorial_null(template: NullConfig, condition: FactorialCondition) -> NullConfig:
    return NullConfig(
        mappability=replace(template.mappability, enabled=condition.mappability),
        gc=replace(template.gc, enabled=condition.gc),
        chromatin=replace(template.chromatin, enabled=condition.chromatin),
    )


# run all eight combinations
def run_factorial(request: AnalysisRequest) -> FactorialResult:
    runs = []

    for condition in FACTORIAL_DESIGN:
        null_config = build_factorial_null(request.null, condition)
        request = replace(request, null=null_config) # preserve other analysis settings

        try:
            result = run_analysis(request)
            runs.append(FactorialRunResult(condition=condition, result=result))
        except (FileNotFoundError, RuntimeError, ValueError) as ex:
            runs.append(FactorialRunResult(condition=condition, error=str(ex)))

    return FactorialResult(runs=runs)
