# src/nullsens/papers/cli.py
# - provides commands to the CLI for using the scripts contained in this folder

# project imports
from nullsens.papers.catalog import refresh_readme
from nullsens.papers.metadata import (
    create_paper, drop_paper, get_existing_paper_ids, get_next_paper_id, load_paper, save_paper
)

# exposed fields for the add/edit cmds
FIELDS = [
    ("title", "Title"),
    ("owner", "Owner"),
    ("doi", "DOI"),
    ("pmid", "PMID"),
    ("data_url", "Data URL"),
    ("code_url", "Code URL"),
    ("reproduction_status", "Reproduction status"),
]

# - command handlers -

# create a new paper record, prompt for metadata, then refresh the catalog
def add_command(args):
    try:
        paper_id = args.paper_id or get_next_paper_id()
        paper = create_paper(paper_id)
    except (ValueError, FileExistsError, RuntimeError) as ex:
        print(f"Error: {ex}")
        return

    print(f"Adding {paper['paper_id']}")
    print("Leave a field blank if it isn't known yet.\n")

    for key, label in FIELDS:
        paper[key] = input(f"{label}: ").strip()

    save_paper(paper["paper_id"], paper)
    refresh_readme()

    print(f"\n{paper['paper_id']} added successfully.")


# update an existing record and then refresh the catalog
def edit_command(args):
    try:
        paper = load_paper(args.paper_id)
    except (ValueError, FileNotFoundError) as ex:
        print(f"Error: {ex}")
        return

    print(f"Editing {paper['paper_id']}")
    print("Press Enter to keep the current value.\n")

    for key, label in FIELDS:
        current = paper.get(key, "")
        value = input(f"{label} [{current}]: ").strip()
        paper[key] = value if value else current

    save_paper(paper["paper_id"], paper)
    refresh_readme()

    print(f"\n{paper['paper_id']} updated successfully.")


# remove a paper record and then refresh the catalog
def delete_command(args):
    try:
        paper_id = drop_paper(args.paper_id)
    except (ValueError, FileNotFoundError) as ex:
        print(f"Error: {ex}")
        return

    refresh_readme()

    print(f"{paper_id} deleted successfully.")


# display all currently assigned paper IDs
def list_command(args):
    paper_ids = get_existing_paper_ids()

    if not paper_ids:
        print("No papers currently exist.")
        return

    print("Existing papers:")
    for paper_id in paper_ids:
        print(f"- {paper_id}")


# regenerate papers/README.md from the YAML metadata files
def refresh_command(args):
    refresh_readme()
    print("Paper README refreshed.")

# - CLI registration -

# register "nullsens papers" and add its subcommands
def register_paper_commands(subparsers):
    papers_parser = subparsers.add_parser("papers", help="Manage project paper metadata")
    paper_commands = papers_parser.add_subparsers(dest="papers_command", required=True)

    # add
    add_parser = paper_commands.add_parser("add", help="Add a new paper")
    add_parser.add_argument("paper_id", nargs="?", help="Optional paper ID, e.g. P01")
    add_parser.set_defaults(func=add_command)

    # edit
    edit_parser = paper_commands.add_parser("edit", help="Edit an existing paper")
    edit_parser.add_argument("paper_id", help="Paper ID, e.g. P01")
    edit_parser.set_defaults(func=edit_command)

    # delete
    delete_parser = paper_commands.add_parser("delete", help="Delete an existing paper")
    delete_parser.add_argument("paper_id", help="Paper ID, e.g. P01")
    delete_parser.set_defaults(func=delete_command)

    # list
    list_parser = paper_commands.add_parser("list", help="List existing paper IDs")
    list_parser.set_defaults(func=list_command)

    # refresh
    refresh_parser = paper_commands.add_parser("refresh", help="Regenerate papers/README.md")
    refresh_parser.set_defaults(func=refresh_command)
