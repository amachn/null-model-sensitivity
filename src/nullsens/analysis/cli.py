# src/nullsens/analysis/cli.py
# - provides commands to the CLI interface for enrichment and factorial analysis

# libraries
from pathlib import Path

# project imports
from nullsens.analysis.factorial import run_factorial
from nullsens.analysis.models import (
    AnalysisRequest,
    ChromatinConfig,
    GCConfig,
    MappabilityConfig,
    NullConfig,
    ReferenceFiles,
)
from nullsens.analysis.reporting import (
    print_analysis_result,
    print_factorial_progress,
    print_factorial_result,
    write_analysis_outputs,
    write_factorial_outputs
)
from nullsens.analysis.runner import run_analysis

# - helpers -

# adds all arguments shared by 'analyze' and 'factorial' commands (prevent duplication)
def add_analysis_arguments(parser, enable_flags: bool) -> None:
    # required analysis inputs
    parser.add_argument("regions_a", type=Path, help="Path to genomic region set A")
    parser.add_argument("regions_b", type=Path, help="Path to genomic region set B")
    parser.add_argument("--genome", required=True, help="Reference genome, e.g. hg38")

    # general analysis options
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use (default: 1)")
    parser.add_argument("--permutations", type=int, default=1000, help="Number of permutations (default: 1000)")
    parser.add_argument("--max-attempts", type=int, default=10000, help="Maximum placement attempts (default: 10000)")
    parser.add_argument("--seed", type=int, help="Optional random seed")
    parser.add_argument("--diagnostics", action="store_true", help="Print backend timing + parallel-exec diagnostics")
    parser.add_argument("--output-prefix", type=Path, help="Write JSON/CSV results using this output path prefix")

    # mappability masking
    if enable_flags:
        parser.add_argument("--mappability", action="store_true", help="Enable mappability masking")

    parser.add_argument("--mappability-reference", type=Path, help="Path to the mappability reference file")
    parser.add_argument("--mappability-threshold", type=float, default=MappabilityConfig().threshold,
                        help=f"Optional mappability score threshold (default: no threshold)")
    parser.add_argument("--mappability-inclusion-workspace", action="store_true",
                        help="Treat the mappability workspace as allowed regions rather than excluded")

    # GC matching
    if enable_flags:
        parser.add_argument("--gc", action="store_true", help="Enable GC-content matching")

    parser.add_argument("--gc-tolerance", type=float, default=GCConfig().tolerance,
                        help=f"Maximum allowed GC-content difference (default: {GCConfig().tolerance})")
    parser.add_argument("--gc-flank", type=int, default=GCConfig().flank,
                        help=f"Number of flanking bases included in calculation (default: {GCConfig().flank})")
    parser.add_argument("--gc-max-n", type=float, default=GCConfig().max_n,
                        help=f"Maximum allowed fraction of N bases (default: {GCConfig().max_n})")

    # chromatin context matching
    if enable_flags:
        parser.add_argument("--chromatin", action="store_true", help="Enable chromatin-context matching")

    parser.add_argument("--chromatin-reference", type=Path, help="Path to the chromatin-state reference file")
    parser.add_argument("--chromatin-padding", type=int, default=ChromatinConfig().padding,
                        help="Padding added around chromatin-state regions "
                            f"(default: {ChromatinConfig().padding})")
    parser.add_argument("--chromatin-label-col", help="Metadata column containing chromatin-state labels")
    parser.add_argument("--chromatin-priority", nargs="+",
                        help="Priority order for overlapping chromatin states, e.g. PLS pELS dELS")
    parser.add_argument("--chromatin-strict", action="store_true",
                        help="Require the entire randomized interval, rather than "
                            "only its midpoint, to remain within the chromatin state")


# builds an AnalysisRequest from parsed CLI arguments
def build_analysis_request(args) -> AnalysisRequest:
    null_config = NullConfig(
        mappability=MappabilityConfig(
            enabled=getattr(args, "mappability", False), # factorial ignores enable/disable flags
            threshold=args.mappability_threshold,
            inclusion_workspace=args.mappability_inclusion_workspace,
        ),
        gc=GCConfig(
            enabled=getattr(args, "gc", False),
            tolerance=args.gc_tolerance,
            flank=args.gc_flank,
            max_n=args.gc_max_n,
        ),
        chromatin=ChromatinConfig(
            enabled=getattr(args, "chromatin", False),
            padding=args.chromatin_padding,
            label_col=args.chromatin_label_col,
            priority=args.chromatin_priority,
            strict=args.chromatin_strict,
        ),
    )

    references = ReferenceFiles(
        mappability=args.mappability_reference,
        chromatin=args.chromatin_reference,
    )

    return AnalysisRequest(
        regions_a=args.regions_a,
        regions_b=args.regions_b,
        genome=args.genome,
        cores=args.cores,
        permutations=args.permutations,
        max_attempts=args.max_attempts,
        seed=args.seed,
        diagnostics=args.diagnostics,
        null=null_config,
        references=references,
    )

# - command handlers -

# run a genomic enrichment analysis using the configured null model
def analyze_command(args):
    request = build_analysis_request(args)

    try:
        print("Running analysis...")
        result = run_analysis(request)
    except (FileNotFoundError, RuntimeError, ValueError) as ex:
        print(f"Error: {ex}")
        return

    print_analysis_result(result)

    if args.output_prefix is not None:
        try:
            json_path, csv_path = write_analysis_outputs(args.output_prefix, request, result)
        except OSError as ex:
            print(f"\nError writing results: {ex}")
            return

        print(f"\nJSON result: {json_path}")
        print(f"CSV result: {csv_path}")


# run a 2^3 factorial enrichment analysis (requires map/chrom refs, runs all null models + combos)
def factorial_command(args):
    request = build_analysis_request(args)

    try:
        result = run_factorial(request, progress=print_factorial_progress)
    except (FileNotFoundError, RuntimeError, ValueError) as ex:
        print(f"Error: {ex}")
        return

    print_factorial_result(result)

    if args.output_prefix is not None:
        try:
            json_path, csv_path = write_factorial_outputs(args.output_prefix, request, result)
        except OSError as ex:
            print(f"\nError writing results: {ex}")
            return

        print(f"\nJSON result: {json_path}")
        print(f"CSV result: {csv_path}")

# - CLI registration -
def register_analysis_commands(subparsers):
    # register "nullsens analyze"
    analyze_parser = subparsers.add_parser("analyze", aliases=["analyse"], help="Run a genomic enrichment analysis")
    add_analysis_arguments(analyze_parser, enable_flags=True)
    analyze_parser.set_defaults(func=analyze_command)

    # register "nullsens factorial"
    factorial_parser = subparsers.add_parser("factorial", help="Run the full factorial combination of null models")
    add_analysis_arguments(factorial_parser, enable_flags=False)
    factorial_parser.set_defaults(func=factorial_command)
