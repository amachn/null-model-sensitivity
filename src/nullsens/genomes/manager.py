# src/nullsens/genomes/manager.py
# - management of R genome (BSgenome) packages required for analyses to be run

# libraries
import subprocess

# project imports
from nullsens.backends.runtime import find_rscript
from nullsens.paths import PACKAGE_PATHS

RSCRIPT = find_rscript()

# - genome queries -

# query installed or available BSgenome packages for a genome
def _query_genomes(genome: str, installed: bool) -> list[str]:
    function = ("installed.genomes" if installed else "available.genomes")
    expression = f"""
        genome <- commandArgs(trailingOnly = TRUE)[1];
        genomes <- BSgenome::{function}(splitNameParts = TRUE);
        
        if (startsWith(genome, "BSgenome.")) {{
            matches <- genomes$pkgname[genomes$pkgname == genome]
        }} else {{
            matches <- genomes$pkgname[genomes$genome == genome & !genomes$masked]
        }};
        
        cat(as.character(matches), sep = "\\n");
    """.replace("\n", "")

    completed = subprocess.run(
        [RSCRIPT, "-e", expression, genome], cwd=PACKAGE_PATHS.r, capture_output=True, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Failed to query BSgenome packages:\n{completed.stderr}")

    return [line for line in completed.stdout.splitlines() if line.strip()]


# resolve a genome assembly or package name to an available BSgenome package
def get_genome_package(genome: str) -> str:
    matches = _query_genomes(genome, installed=False)

    if not matches:
        raise ValueError(f"No available BSgenome package found for genome '{genome}'.")

    if len(matches) > 1:
        packages = "\n".join(f" - {package}" for package in matches)
        raise ValueError(
            f"Genome '{genome}' is ambiguous.\n"
            f"Matching packages:\n{packages}\n"
            "Specify the full BSgenome package name."
        )

    return matches[0]


# check whether a genome is installed in the nullsens R environment
def is_genome_installed(genome: str) -> bool:
    matches = _query_genomes(genome, installed=True)

    if len(matches) > 1:
        packages = "\n".join(f" - {package}" for package in matches)
        raise ValueError(
            f"Genome '{genome}' is ambiguous.\n"
            f"Installed matching packages:\n{packages}"
        )

    return len(matches) == 1


# list BSgenome packages installed in the nullsens R environment
def get_installed_genomes() -> list[tuple[str, str]]:
    expression = """
        genomes <- BSgenome::installed.genomes(splitNameParts = TRUE);
        if (nrow(genomes) > 0) {
            cat(paste(genomes$genome, genomes$pkgname, sep = "\\t"), sep = "\\n")  
        };
    """.replace("\n", "")

    completed = subprocess.run(
        [RSCRIPT, "-e", expression], cwd=PACKAGE_PATHS.r, capture_output=True, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Failed to query installed genomes:\n{completed.stderr}")

    genomes = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue

        genome, package = line.split("\t", 1)
        genomes.append((genome, package))

    return genomes

# - genome management -

# install a BSgenome package into the nullsens R environment
def install_genome(genome: str) -> str:
    package = get_genome_package(genome)
    if is_genome_installed(package):
        return package

    expression = f"""
        package <- commandArgs(trailingOnly = TRUE)[1];
        renv::install(paste0("bioc::", package), prompt = FALSE);
    """.replace("\n", "")

    completed = subprocess.run(
        [RSCRIPT, "-e", expression, package], cwd=PACKAGE_PATHS.r, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Failed to install genome '{genome}'.")

    return package


# remove a BSgenome package from the nullsens R environment
def remove_genome(genome: str) -> str:
    matches = _query_genomes(genome, installed=True)

    if not matches:
        raise ValueError(f"Genome '{genome}' is not installed.")

    if len(matches) > 1:
        packages = "\n".join(f" - {package}" for package in matches)
        raise ValueError(
            f"Genome '{genome}' is ambiguous.\n"
            f"Installed matching packages:\n{packages}\n"
            "Specify the full BSgenome package name."
        )

    package = matches[0]
    expression = """
        package <- commandArgs(trailingOnly = TRUE)[1];
        renv::remove(package, prompt = FALSE);
    """.replace("\n", "")

    completed = subprocess.run(
        [RSCRIPT, "-e", expression, package], cwd=PACKAGE_PATHS.r, text=True
    )

    if completed.returncode != 0:
        raise RuntimeError(f"Failed to remove genome '{genome}'.")

    return package
