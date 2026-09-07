# src/nullsens/cli.py
# - primary command line interface bridge for functions found in the source code

# libraries
import argparse

# project imports
from nullsens.papers.cli import register_paper_commands


# CLI entrypoint
def main():
    # top-level CLI parser for the "nullsens" command
    parser = argparse.ArgumentParser(prog="nullsens", description="Null Model Sensitivity project tools")

    # register first-level subcommands, e.g. "nullsens papers"
    subparsers = parser.add_subparsers(dest="command", required=True)

    register_paper_commands(subparsers)

    # parse CLI input and dispatch to the selected command
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()