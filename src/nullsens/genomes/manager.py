# src/nullsens/genomes/manager.py
# - management of R genome (BSgenome) packages required for analyses to be run

# libraries
import subprocess

# project imports
from nullsens.backends.r_runtime import find_rscript
from nullsens.paths import PACKAGE_PATHS


SUPPORTED_GENOMES = {
    "hg38": "BSgenome.Hsapiens.UCSC.hg38",
    "mm10": "BSgenome.Mmusculus.UCSC.mm10",
}


def get_genome_package(genome: str) -> str:
    try:
        return SUPPORTED_GENOMES[genome]
    except KeyError:
        supported = ", ".join(SUPPORTED_GENOMES)
        raise ValueError(f"Unsupported genome '{genome}'.\nSupported genomes: {supported}")


def is_genome_installed(genome: str) -> bool:
    package = get_genome_package(genome)
    rscript = find_rscript()

    expression = f"quit(status = if(system.file(package = '{package}') == '') 1 else 0)"
    completed = subprocess.run([rscript, "-e", expression], cwd=PACKAGE_PATHS.r, capture_output=True, text=True)

    return completed.returncode == 0
