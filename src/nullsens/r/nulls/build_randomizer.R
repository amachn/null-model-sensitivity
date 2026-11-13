# ~~~ src/nullsens/r/nulls/build_randomizer.R ~~~
# - generates a randomization function appropriate for the provided null combo

# evaluate every enabled validator
validate_candidates <- function(original, candidate, validators) {
  valid <- rep(TRUE, length(candidate))
  for (validator in validators) {
    result <- validator(original, candidate)
    if (!is.logical(result) || length(result) != length(candidate)) {
      stop("A validator must return one TRUE/FALSE per candidate region.")
    }
    valid <- valid & result
  }
  valid
}

build_randomizer <- function(
  config, refs, regions, genome, max_attempts = 100L
) {
  mask <- NA
  validators <- list()
  
  # mappability contributes a forbidden-region mask
  if (config$mappability_masked) {
    if (is.null(refs$mappability)) {
      stop("Mappability masking is enabled but no reference was provided.")
    }
    
    mask <- prepare_mappability(
      refs$mappability, genome, 
      config$mappability$threshold, config$mappability$invert
    )
  }
  
  # GC and chromatin modules will add validators here
  if (config$gc_matched) {
    stop("GC matching is enabled but hasn't been integrated yet.")
  }
  
  if (config$chromatin_matched) {
    stop("Chromatin matching is enabled but hasn't been integrated yet.")
  }
  
  function(A, genome, ...) {
    # TODO: will regioneR pass a mask itself, and if so do we need to combine?
    randomized <- regioneR::randomizeRegions(
      A = a, genome = genome, mask = mask, ...
    )
    
    # static masks require no resampling
    if (length(validators) == 0L) {
      return(randomized)
    }
    
    valid <- validate_candidates(A, randomized, validators)
    attempts <- rep(1L, length(A))
    
    # resample only regions that fail the enabled constraints
    while (any(!valid)) {
      failed <- which(!valid)
      
      attempts[failed] <- attempts[failed] + 1L
      exceeded <- failed[attempts[failed] > max_attempts]
      
      if (length(exceeded) > 0L) {
        stop(
          "Unable to place ", length(exceeded), " region(s) after ",
          max_attempts, " attempts under the requested null-model constraints."
        )
      }
      
      replacements <- regioneR::randomizeRegions(
        A = A[failed],
        genome = genome,
        mask = mask,
        ...
      )
      
      replacement_valid <- validate_candidates(
        A[failed],
        replacements,
        validators
      )
      
      accepted <- failed[replacement_valid]

      # TODO: need to finalize adding the valid replacements and updating the valid vector
    }
    
    randomized
  }
}

#Anugraka - Wrote a version of the build randomizer function to test map + gc, not sure if this is what we wanted
'''
build_randomizer <- function(
  config, refs, regions, genome, max_attempts = 100L
) {
  mask <- NULL
  map_ctx <- NULL
  validators <- list()
  diag <- new.env(parent = emptyenv())
  diag$attempts <- integer(0)

  if (isTRUE(config$mappability_masked)) {
    if (is.null(refs$mappability)) {
      stop("Mappability masking is enabled but no reference was provided.")
    }

    map_ctx <- prepare_mappability(
      refs$mappability, genome,
      config$mappability$threshold, config$mappability$invert
    )
    mask <- map_ctx$mask
  }

  if (isTRUE(config$gc_matched)) {
    gc_ctx <- prepare_gc(
      genome,
      tolerance = config$gc$tolerance,
      flank     = config$gc$flank,
      max_n     = config$gc$max_n
    )
    validators <- c(validators, list(function(original, candidate) {
      validate_gc(gc_ctx, original, candidate)
    }))
  }

  if (isTRUE(config$chromatin_matched)) {
    stop("Chromatin matching is enabled but hasn't been integrated yet.")
  }

  function(A, genome, ...) {
    randomized <- regioneR::randomizeRegions(
      A = A, genome = genome, mask = mask, ...
    )

    if (length(validators) == 0L) {
      diag$attempts <- rep(1L, length(A))
      return(randomized)
    }

    valid <- validate_candidates(A, randomized, validators)
    attempts <- rep(1L, length(A))

    while (any(!valid)) {
      failed <- which(!valid)

      attempts[failed] <- attempts[failed] + 1L
      exceeded <- failed[attempts[failed] > max_attempts]

      if (length(exceeded) > 0L) {
        diag$attempts <- attempts
        stop(
          "Unable to place ", length(exceeded), " region(s) after ",
          max_attempts, " attempts under the requested null-model constraints."
        )
      }

      replacements <- regioneR::randomizeRegions(
        A = A[failed],
        genome = genome,
        mask = mask,
        ...
      )

      replacement_valid <- validate_candidates(
        A[failed],
        replacements,
        validators
      )

      accepted <- failed[replacement_valid]
      randomized[accepted] <- replacements[replacement_valid]
      valid[accepted] <- TRUE
    }

    diag$attempts <- attempts
    randomized
  }
}
'''

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
