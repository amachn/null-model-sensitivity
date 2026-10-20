# ~~~ src/nullsens/r/nulls/length.R ~~~
# - implementation of the length-matching null model

# randomize positions while preserving the exact width of each original interval
randomize_exact_length <- function(
    A, genome, mask = NULL, ...
) {
  regioneR::randomizeRegions(A = a, genome = genome, mask = mask)
}

# placeholder
randomize_length_distribution <- function(
    A, genome, mask = NULL, ...
) {
  stop("length-distribution matching isn't implemented yet.")
}

# return the randomization function according to the length mode
get_length_randomizer <- function(mode) {
  if (mode == "exact") { return(randomize_exact_length) }
  if (mode == "distribution") { return(randomize_length_distribution) }
  stop("unknown length mode: ", mode)
}