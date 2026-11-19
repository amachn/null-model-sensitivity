# src/nullsens/analysis/cli.py
# - provides commands to the CLI interface for analysis

# libraries
from pathlib import Path

# project imports
from nullsens.analysis.models import (
    AnalysisRequest,
    ChromatinConfig,
    GCConfig,
    MappabilityConfig,
    NullConfig,
    ReferenceFiles,
)
from nullsens.analysis.runner import run_analysis

# - command handlers -

# run a genomic enrichment analysis using the configured null model
def analyze_command(args):
    null_config = NullConfig(
        chromatin=ChromatinConfig(
            enabled=args.chromatin,
            padding=args.chromatin_padding,
            label_col=args.chromatin_label_col,
            priority=args.chromatin_priority,
            strict=args.chromatin_strict,
        ),
        gc=GCConfig(
            enabled=args.gc,
            tolerance=args.gc_tolerance,
            flank=args.gc_flank,
            max_n=args.gc_max_n,
        ),
        mappability=MappabilityConfig(
            enabled=args.mappability,
            threshold=args.mappability_threshold,
            inclusion_workspace=args.mappability_inclusion_workspace,
        ),
    )

    references = ReferenceFiles(
        chromatin=args.chromatin_reference,
        mappability=args.mappability_reference,
    )

    request = AnalysisRequest(
        regions_a=args.regions_a,
        regions_b=args.regions_b,
        genome=args.genome,
        permutations=args.permutations,
        max_attempts=args.max_attempts,
        seed=args.seed,
        null=null_config,
        references=references,
    )

    try:
        result = run_analysis(request)
    except (FileNotFoundError, RuntimeError, ValueError) as ex:
        print(f"Error: {ex}")
        return

    print("Analysis complete.")
    print(f"Observed overlap: {result.observed:g}")
    print(f"Null mean: {result.null_mean:.4f}")
    print(f"Null SD: {result.null_sd:.4f}")
    print(f"Z-score: {result.z_score:.4f}")
    print(f"P-value: {result.p_value:.6g}")

# - CLI registration -

# register "nullsens analyze"
def register_analysis_commands(subparsers):
    analyze_parser = subparsers.add_parser("analyze", aliases=["analyse"], help="Run a genomic enrichment analysis")

    # required analysis inputs
    analyze_parser.add_argument("regions_a", type=Path, help="Path to genomic region set A")
    analyze_parser.add_argument("regions_b", type=Path, help="Path to genomic region set B")
    analyze_parser.add_argument("--genome", required=True, help="Reference genome, e.g. hg38")

    # general analysis options
    analyze_parser.add_argument("--permutations", type=int, default=1000,
                                help="Number of permutations (default: 1000)")
    analyze_parser.add_argument("--max-attempts", type=int, default=100,
                                help="Maximum placement attempts per region (default: 100)")
    analyze_parser.add_argument("--seed", type=int, help="Optional random seed")

    # mappability masking
    analyze_parser.add_argument("--mappability", action="store_true", help="Enable mappability masking")
    analyze_parser.add_argument("--mappability-reference", type=Path, help="Path to the mappability reference file")
    analyze_parser.add_argument("--mappability-threshold", type=float, default=MappabilityConfig().threshold,
                                help=f"Optional mappability score threshold (default: {MappabilityConfig().threshold})")
    analyze_parser.add_argument("--mappability-inclusion-workspace", action="store_true",
                                help="Treat the mappability workspace as allowed regions rather than excluded")

    # GC matching
    analyze_parser.add_argument("--gc", action="store_true", help="Enable GC-content matching")
    analyze_parser.add_argument("--gc-tolerance", type=float, default=GCConfig().tolerance,
                                help=f"Maximum allowed GC-content difference (default: {GCConfig().tolerance})")
    analyze_parser.add_argument("--gc-flank", type=int, default=GCConfig().flank,
                                help=f"Number of flanking bases included in calculation (default: {GCConfig().flank})")
    analyze_parser.add_argument("--gc-max-n", type=float, default=GCConfig().max_n,
                                help=f"Maximum allowed fraction of N bases (default: {GCConfig().max_n})")

    # chromatin-context matching
    analyze_parser.add_argument("--chromatin", action="store_true", help="Enable chromatin-context matching")
    analyze_parser.add_argument("--chromatin-reference", type=Path, help="Path to the chromatin-state reference file")
    analyze_parser.add_argument("--chromatin-padding", type=int, default=ChromatinConfig().padding,
                                help="Padding added around chromatin-state regions "
                                     f"(default: {ChromatinConfig().padding})")
    analyze_parser.add_argument("--chromatin-label-col", help="Metadata column containing chromatin-state labels")
    analyze_parser.add_argument("--chromatin-priority", nargs="+",
                                help="Priority order for overlapping chromatin states, e.g. PLS pELS dELS")
    analyze_parser.add_argument("--chromatin-strict", action="store_true",
                                help="Require the entire randomized interval, rather than "
                                     "only its midpoint, to remain within the chromatin state")

    analyze_parser.set_defaults(func=analyze_command)