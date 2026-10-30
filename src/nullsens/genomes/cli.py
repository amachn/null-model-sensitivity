# src/nullsens/genomes/cli.py
# - provides commands to the CLI interface for managing reference genomes

# project imports
from nullsens.genomes.manager import get_installed_genomes, install_genome, remove_genome

# - command handlers -

# list installed reference genomes
def list_command(args):
    try:
        genomes = get_installed_genomes()
    except RuntimeError as ex:
        print(f"Error: {ex}")
        return

    if not genomes:
        print("No reference genomes installed.")
        return

    for genome, package in genomes:
        print(f"{genome}: {package}")


# install a reference genome
def install_command(args):
    try:
        package = install_genome(args.genome)
    except (ValueError, RuntimeError) as ex:
        print(f"Error: {ex}")
        return

    print(f"Genome '{args.genome}' installed ({package}).")


# remove a reference genome
def remove_command(args):
    try:
        package = remove_genome(args.genome)
    except (ValueError, RuntimeError) as ex:
        print(f"Error: {ex}")
        return

    print(f"Genome '{args.genome}' removed ({package}).")

# - CLI registration -
def register_genome_commands(subparsers):
    genomes_parser = subparsers.add_parser("genomes", help="Manage reference genomes")
    genome_commands = genomes_parser.add_subparsers(dest="genome_command", required=True)

    # install
    install_parser = genome_commands.add_parser("install", help="Install a reference genome")
    install_parser.add_argument("genome", help="Genome assembly or full BSgenome package name")
    install_parser.set_defaults(func=install_command)

    # remove
    remove_parser = genome_commands.add_parser("remove", help="Remove an installed reference genome")
    remove_parser.add_argument("genome", help="Genome assembly or full BSgenome package name")
    remove_parser.set_defaults(func=remove_command)

    # list
    list_parser = genome_commands.add_parser("list", help="List installed reference genomes")
    list_parser.set_defaults(func=list_command)