# ~~~ src/nullsens/r/nulls/build_randomizer.R ~~~
# - generates a randomization function appropriate for the provided null combo
# * INTEGRATION OF THE NULLS HAPPENS HERE
source("length.R")

build_randomizer <- function(config, refs, regions, genome) {
  # current thoughts:
  # 1. establish length behavior (operates differently from other nulls)
  # 2. load/apply requested constraints (gc, chrom, mappability, etc)
  # 3. combine into a regioneR-compatible randomization function
  #    - currently unsure of how we'll do this
  # 4. return that function
  get_length_randomizer(config$length_mode)
}