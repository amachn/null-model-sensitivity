# ~~~ src/nullsens/r/nulls/build_randomizer.R ~~~
# - generates a randomization function appropriate for the provided null combo
# * INTEGRATION OF THE NULLS HAPPENS HERE
build_randomizer <- function(config, refs, regions, genome) {
  # core responsibilities
  # - at build time:
  # prepare GC context (if enabled)
  # prepare mappability context (if enabled)
  # prepare chromatin context (if enabled)
  # return randomization function (as a function closure, which does the "permutation time" tasks listed below)
  # - at permutation time (for each perm):
  # unresolved = every region
  # while unresolved regions exist: (needs a max attempt limit)
  #   randomize unresolved regions (preserving exact widths/lengths [baseline model, always on])
  #   test proposals against: mappability, GC, chromatin (based on config)
  #   accept valid candidates
  #   unresolved = failed proposals
  # return complete randomized GRanges

  # notes:
  # - need to preserve ordering/numbering of the sequences when shuffling, i.e., original[1] == randomized[1]
  #   - necessary for chromatin/GC matching as they are dependent on the original interval
  # - nulls maybe follow a two-function format of:
  #   - prepare_<null> -> loads/preprocesses reference data, one-time only at build time
  #   - validate_<null> -> takes original and candidate regions and returns a logical vector like:
  #     [TRUE TRUE FALSE TRUE FALSE ...]
  #     indicating which candidates satisfy the constraints applied by the enabled nulls
}