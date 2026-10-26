# src/nullsens/analysis/cli.py
# - provides commands to the CLI interface for analysis

# libraries
from pathlib import Path

# project imports
from nullsens.analysis.models import AnalysisRequest
from nullsens.analysis.runner import run_analysis

# - command handlers -

# run a genomic enrichment analysis using the configured null model
def analyze_command(args):
    request = AnalysisRequest(
        regions_a=args.regions_a,
        regions_b=args.regions_b,
        genome=args.genome,
        permutations=args.permutations,
        seed=args.seed
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

    # arguments
    analyze_parser.add_argument("regions_a", type=Path, help="Path to genomic region set A")
    analyze_parser.add_argument("regions_b", type=Path, help="Path to genomic region set B")
    analyze_parser.add_argument("--genome", required=True, help="Reference genome, e.g. hg38")
    analyze_parser.add_argument("--permutations", type=int, default=1000, help="Number of permutations (default: 1000)")
    analyze_parser.add_argument("--seed", type=int, help="Optional random seed")

    # TODO: expand arguments with null flags once nulls are integrated

    analyze_parser.set_defaults(func=analyze_command)