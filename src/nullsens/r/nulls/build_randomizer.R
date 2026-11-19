# ~~~ src/nullsens/r/nulls/build_randomizer.R ~~~
# - builds one regioneR-compatible randomization function from the
#   requested null constraints
#
# - at build time:
#   - prepare mappability mask, if enabled
#   - prepare GC context, if enabled
#   - prepare chromatin context, if enabled
#   - construct validator functions
# - at permutation time:
#   - generate exact-width randomized regions
#   - enforce static workspace masks during randomization
#   - validate GC/chromatin constraints
#   - resample only failed regions until every region passes

# runs every enabled per-region validator
# - each validator must return exactly one T/F value for each candidate region:
#   - TRUE = candidate satisfies the constraint
#   - FALSE = candidate must be resampled
# - a candidate is accepted only when it passes every enabled validator
validate_candidates <- function(original, candidate, validators) {
  valid <- rep(TRUE, length(candidate))
  
  for (name in names(validators)) {
    result <- validators[[name]](original, candidate)
    
    if (!is.logical(result) || length(result) != length(candidate)) {
      stop(
        "Null validator '", name, "' must return one TRUE/FALSE value ",
        "per candidate region."
      )
    }
    
    # an indeterminate result can't be accepted
    result[is.na(result)] <- FALSE
    valid <- valid & result
  }
  
  valid
}

# verify invariants required of every nullsens randomization
# - the randomized set must:
#   - contain exactly one candidate per original region
#   - preserve the exact width of each corresponding original region
# - maintaining index order is important because GC/chromatin matching
#   compare candidate[i] against original[i]
validate_randomization <- function(original, randomized) {
  if (length(original) != length(randomized)) {
    stop("Randomization changed the number of regions.")
  }
  
  if (
    !all(GenomicRanges::width(original) == GenomicRanges::width(randomized))
  ) {
    stop("Randomization failed to preserve exact region widths.")
  }
  
  invisible(TRUE)
}

# build one randomization function for one null-model config
# - `regions` is kept because it is the original region set
#   - future impls may use it to precompute original-region properties once
#     rather than repeatedly
build_randomizer <- function(
  config, refs, regions, genome, max_attempts = 100L
) {
  if (max_attempts < 1L) { stop("max_attempts must be an integer >= 1.") }
  
  # ----
  # 1 -> prepare context, static mask and candidate validators
  # ----
  # these context objects need to exist outside of their if statements
  mappability_context <- NULL
  chromatin_context <- NULL
  
  # NA means nullsens explicitly has no baseline forbidden mask
  proposal_mask <- NA
  
  # GC and chrom validators will be stored here
  validators <- list()
  
  # ----
  # 2 -> prepare mappability
  # ----
  # mappability is primarily a static constraint rather than a validator.
  # its forbidden regions can thus be handed directly to regioneR.
  if (isTRUE(config$mappability$enabled)) {
    if (is.null(refs$mappability)) {
      stop("Mappability masking is enabled but no reference was provided.")
    }
    
    mappability_context <- prepare_mappability(
      mappability_bed     = refs$mappability,
      genome              = genome,
      threshold           = config$mappability$threshold,
      inclusion_workspace = config$mappability$inclusion_workspace
    )
    
    # normally, mappability is enforced directly by randomizeRegions()
    proposal_mask <- mappability_context$mask
  }
  
  # ----
  # 3 -> prepare GC matching
  # ----
  # gc is region-specific: 
  # candidate[i] must have approx. the same GC fraction as original[i].
  # because of this, it becomes a validator used during rejection/resampling.
  if (isTRUE(config$gc$enabled)) {
    gc_context <- prepare_gc(
      genome,
      tolerance = config$gc$tolerance,
      flank     = config$gc$flank,
      max_n     = config$gc$max_n
    )
    
    validators$gc <- function(original, candidate) {
      validate_gc(
        context   = gc_context, 
        original  = original,
        candidate = candidate
      )
    }
  }
  
  # ----
  # 4 -> prepare chromatin-context matching
  # ----
  # like gc, chromatin depends on the corresponding original interval, so it
  # is also represented by a validator rather than a static mask.
  if (isTRUE(config$chromatin$enabled)) {
    if (is.null(refs$chromatin)) {
      stop("Chromatin matching is enabled but no reference was provided.")
    }
    
    chromatin_context <- prepare_chromatin(
      refs$chromatin, genome,
      padding   = config$chromatin$padding,
      label_col = config$chromatin$label_col,
      priority  = config$chromatin$priority
    )
    
    # keep validation even though the proposer is already chromatin-aware.
    validators$chromatin <- function(original, candidate) {
      validate_chromatin(
        context   = chromatin_context,
        original  = original,
        candidate = candidate,
        strict    = isTRUE(config$chromatin$strict)
      )
    }
    
    # propose_chromatin() doesn't understand regioneR masks - thus, when
    # mappability + chromatin are both enabled, mappability must move from
    # proposal-time enforcement to post-proposal validation
    if (!is.null(mappability_context)) {
      validators$mappability <- function(original, candidate) {
        validate_mappability(
          context   = mappability_context,
          original  = original,
          candidate = candidate
        )
      }
      
      # the chromatin proposer will be used instead of randomizeRegions(),
      # so this static mask is no longer used
      proposal_mask <- NA
    }
  }
  
  # ----
  # 5 -> define the randomization proposer
  # ----
  propose <- function(A, genome, ...) {
    # chromatin has its own state-aware proposal mechanism
    if (!is.null(chromatin_context)) {
      return(propose_chromatin(chromatin_context, A))
    }
    
    # otherwise, use regioneR's normal exact-width randomization
    args <- list(...)
    args$mask <- NULL # ignore any externally supplied mask
    
    do.call(
      regioneR::randomizeRegions,
      c(list(A = A, genome = genome, mask = proposal_mask), args)
    )
  }
  
  # everything above happens ONCE when the analysis is constructed
  # ----
  # 6 -> return the regioneR-compatible randomization closure
  # ----
  # everything below happens repeatedly, once per permutation
  function(A, genome, mask = NA, ...) {
    # -----
    # 6A -> generate one candidate for every original region
    # -----
    randomized <- propose(A = A, genome = genome, ...)
    validate_randomization(A, randomized)
    
    # -----
    # 6B -> if no region-specific constraints are enabled, finish here
    # -----
    if (length(validators) == 0L) {
      return(randomized)
    }
    
    # -----
    # 6C -> check the first set of proposals
    # -----
    valid <- validate_candidates(
      original = A, candidate = randomized, validators = validators
    )
    
    # initial proposal counts as attempt 1
    attempts <- rep(1L, length(A))
    
    # -----
    # 6D -> resample only failed candidates
    # -----
    while(any(!valid)) {
      failed <- which(!valid)
      
      # stop before making another proposal for regions that have already 
      # reached the allowed attempt count
      exceeded <- failed[attempts[failed] >= max_attempts]
      if (length(exceeded) > 0L) {
        stop(
          "Unable to place ", length(exceeded), " region(s) after ",
          max_attempts, " attempts under the requested null-model constraints."
        )
      }
      
      # increment only the regions that need another proposal
      attempts[failed] <- attempts[failed] + 1L
      
      # generate replacements for failed regions only
      replacements <- propose(A = A[failed], genome = genome, ...)
      validate_randomization(A[failed], replacements)
      
      # test the replacement candidates against every enabled constraint
      replacement_valid <- validate_candidates(
        original = A[failed], candidate = replacements, validators = validators
      )
      
      # convert successful replacement positions back into indices within
      # the complete randomized region set
      accepted <- failed[replacement_valid]
      
      # save only successful replacements
      if (length(accepted) > 0L) {
        randomized[accepted] <- replacements[replacement_valid]
        valid[accepted] <- TRUE
      }
    }
    
    # -----
    # 6E -> final invariant check and return
    # -----
    validate_randomization(A, randomized)
    
    randomized
  }
}
