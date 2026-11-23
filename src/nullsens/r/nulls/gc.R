#original gc.R
'''
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Biostrings)
})

prepare_gc <- function(genome, tolerance = 0.02, flank = 0L, max_n = 0.1) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0) {
    stop("GC tolerance must be a single number > 0, got: ", paste(tolerance, collapse = ", "))
  }
  if (!is.numeric(flank) || length(flank) != 1L || flank < 0) {
    stop("GC flank must be a single number >= 0, got: ", paste(flank, collapse = ", "))
  }
  if (!is.numeric(max_n) || length(max_n) != 1L || max_n < 0 || max_n > 1) {
    stop("GC max_n must be between 0 and 1, got: ", paste(max_n, collapse = ", "))
  }

  # GC is computed from the genome sequence, so a BSgenome is required
  sequences <- if (inherits(genome, "BSgenome") || methods::is(genome, "DNAStringSet")) {
    genome
  } else {
    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("GC matching requires the BSgenome package.")
    }
    tryCatch(
      suppressWarnings(
        suppressMessages(
          BSgenome::getBSgenome(genome, masked = FALSE)
        )
      ),
      error = function(err) {
        stop(
          "GC matching could not load a BSgenome for '", genome, "'.\n",
          conditionMessage(err), "\n",
          "Run: nullsens genomes install ", genome
        )
      }
    )
  }

  genome_seqs <- if (inherits(sequences, "BSgenome")) {
    GenomeInfoDb::seqlengths(sequences)
  } else {
    stats::setNames(as.numeric(Biostrings::width(sequences)), names(sequences))
  }

  # chromosome sequences are loaded on first use and cached for later permutations
  cache <- new.env(parent = emptyenv())
  fetch <- function(chr) {
    if (is.null(cache[[chr]])) {
      seq <- sequences[[chr]]
      if (methods::is(seq, "MaskedDNAString")) seq <- Biostrings::unmasked(seq)
      cache[[chr]] <- seq
    }
    cache[[chr]]
  }

  list(
    fetch       = fetch,
    seqnames    = names(genome_seqs),
    source      = if (is.character(genome)) genome else class(genome)[1],
    tolerance   = tolerance,
    flank       = as.integer(flank),
    max_n       = max_n,
    n_sequences = length(genome_seqs)
  )
}

# GC fraction and N fraction of each interval (+/- flank), order preserved
gc_content <- function(context, regions) {
  n <- length(regions)
  gc <- rep(NA_real_, n)
  n_frac <- rep(1, n)

  chroms <- as.character(GenomicRanges::seqnames(regions))
  starts <- GenomicRanges::start(regions) - context$flank
  ends <- GenomicRanges::end(regions) + context$flank

  for (chr in intersect(unique(chroms), context$seqnames)) {
    idx <- which(chroms == chr)
    seq <- context$fetch(chr)
    views <- IRanges::Views(
      seq,
      start = as.integer(pmax(1, starts[idx])),
      end = as.integer(pmin(length(seq), ends[idx]))
    )
    widths <- IRanges::width(views)
    n_count <- Biostrings::letterFrequency(views, "N")[, 1]
    gc_count <- Biostrings::letterFrequency(views, "GC")[, 1]
    called <- widths - n_count

    gc[idx] <- ifelse(called > 0, gc_count / called, NA_real_)
    n_frac[idx] <- n_count / widths
  }

  list(gc = gc, n_frac = n_frac)
}

validate_gc <- function(context, original, candidate) {
  if (length(original) != length(candidate)) {
    stop("validate_gc: original and candidate must have the same length and order.")
  }

  orig <- gc_content(context, original)
  cand <- gc_content(context, candidate)

  valid <- abs(cand$gc - orig$gc) <= context$tolerance & cand$n_frac <= context$max_n

  # an original with no called bases has no GC to match; any low-N candidate is fine
  no_gc <- is.na(orig$gc)
  valid[no_gc] <- cand$n_frac[no_gc] <= context$max_n

  valid[is.na(valid)] <- FALSE
  valid
}
'''

#new gc.R
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Biostrings)
})

# C++ helpers (compiled on first use)
gc_index_cpp <- '
#include <Rcpp.h>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

struct GcIdx { const Rbyte* words; const int* prefix; const int* ns; const int* ne; const int* nc; int nr; int len; };

static std::vector<GcIdx> gc_unpack(const List& idx) {
  std::vector<GcIdx> v(idx.size());
  for (R_xlen_t i = 0; i < idx.size(); ++i) {
    List e = idx[i];
    RawVector w = e["words"]; IntegerVector pf = e["prefix"];
    IntegerVector ns = e["nstart"], ne = e["nend"], nc = e["ncum"];
    GcIdx g; g.words = w.begin(); g.prefix = pf.begin(); g.ns = ns.begin(); g.ne = ne.begin(); g.nc = nc.begin();
    g.nr = ns.size(); g.len = as<int>(e["len"]);
    v[i] = g;
  }
  return v;
}
// G/C count before position x
static inline int gc_before(const GcIdx& g, int x) {
  int w = x >> 6, r = x & 63;
  int c = g.prefix[w];
  if (r) { uint64_t v; std::memcpy(&v, g.words + (size_t)w * 8, 8); c += __builtin_popcountll(v & ((1ULL << r) - 1ULL)); }
  return c;
}
// N count before position x
static inline int n_before(const GcIdx& g, int x) {
  if (g.nr == 0) return 0;
  int j = (int)(std::lower_bound(g.ns, g.ns + g.nr, x) - g.ns);
  if (j == 0) return 0;
  int k = j - 1;
  return g.nc[k] + std::min(x, g.ne[k]) - g.ns[k];
}
// does one window pass the validate_gc test? (target NaN: only the N check applies)
static inline bool window_ok(const GcIdx& c, int start, int width, double target, double tol, double max_n, int flank) {
  int l = std::max(1, start - flank);
  int r = std::min(c.len, start + width - 1 + flank);
  if (l > r) return false;
  double win = (double)(r - l + 1);
  double nct = (double)(n_before(c, r) - n_before(c, l - 1));
  if (!((nct / win) <= max_n)) return false;
  if (std::isnan(target)) return true;
  double called = win - nct;
  if (called <= 0) return false;
  double gcct = (double)(gc_before(c, r) - gc_before(c, l - 1));
  return std::fabs(gcct / called - target) <= tol;
}

// [[Rcpp::export]]
List gc_build_index(RawVector seq) {
  const int L = seq.size();
  const int nw = (L + 63) / 64;
  bool is_gc[256], is_n[256];
  for (int i = 0; i < 256; ++i) { is_gc[i] = false; is_n[i] = false; }
  is_gc[71] = is_gc[103] = is_gc[67] = is_gc[99] = true;
  is_n[78] = is_n[110] = true;
  std::vector<uint64_t> buf(nw, 0ULL);
  std::vector<int> ns, ne;
  int run_start = -1;
  const Rbyte* s = seq.begin();
  for (int i = 0; i < L; ++i) {
    unsigned char b = s[i];
    if (is_gc[b]) buf[i >> 6] |= (1ULL << (i & 63));
    if (is_n[b]) { if (run_start < 0) run_start = i; }
    else if (run_start >= 0) { ns.push_back(run_start); ne.push_back(i); run_start = -1; }
  }
  if (run_start >= 0) { ns.push_back(run_start); ne.push_back(L); }
  IntegerVector prefix(nw + 1);
  int acc = 0;
  for (int w = 0; w < nw; ++w) { prefix[w] = acc; acc += __builtin_popcountll(buf[w]); }
  prefix[nw] = acc;
  RawVector words((R_xlen_t)nw * 8);
  if (nw > 0) std::memcpy(words.begin(), buf.data(), (size_t)nw * 8);
  IntegerVector nstart(ns.begin(), ns.end()), nend(ne.begin(), ne.end()), ncum(ns.size());
  int cum = 0;
  for (size_t k = 0; k < ns.size(); ++k) { ncum[k] = cum; cum += ne[k] - ns[k]; }
  return List::create(_["words"] = words, _["prefix"] = prefix, _["nstart"] = nstart, _["nend"] = nend,
                      _["ncum"] = ncum, _["len"] = L);
}

// G/C and N counts for each window (1-based, inclusive)
// [[Rcpp::export]]
List gc_window_stats(List idx, IntegerVector chr_id, IntegerVector left, IntegerVector right) {
  std::vector<GcIdx> g = gc_unpack(idx);
  const R_xlen_t n = chr_id.size();
  NumericVector gc(n), nn(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const GcIdx& c = g[chr_id[i] - 1];
    int l = left[i] - 1, r = right[i];
    gc[i] = gc_before(c, r) - gc_before(c, l);
    nn[i] = n_before(c, r) - n_before(c, l);
  }
  return List::create(_["gc"] = gc, _["n"] = nn);
}

// samplers only: which windows pass the test?
// [[Rcpp::export]]
LogicalVector gc_window_ok(List idx, IntegerVector chr_id, IntegerVector start, IntegerVector width,
                           NumericVector target, double tol, double max_n, int flank) {
  std::vector<GcIdx> g = gc_unpack(idx);
  const R_xlen_t n = chr_id.size();
  LogicalVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = window_ok(g[chr_id[i] - 1], start[i], width[i], target[i], tol, max_n, flank);
  }
  return out;
}

// samplers only: all passing starts on one chromosome
// [[Rcpp::export]]
List gc_all_compatible_starts(List idx, int chr, int first_start, int last_start, int width,
                              double target, double tol, double max_n, int flank, int cap) {
  std::vector<GcIdx> g = gc_unpack(idx);
  const GcIdx& c = g[chr - 1];
  std::vector<int> out;
  bool overflow = false;
  for (int start = first_start; start <= last_start; ++start) {
    if (window_ok(c, start, width, target, tol, max_n, flank)) {
      if ((int)out.size() >= cap) { overflow = true; break; }
      out.push_back(start);
    }
  }
  return List::create(_["starts"] = IntegerVector(out.begin(), out.end()), _["overflow"] = overflow);
}
'
gc_lib <- new.env(parent = emptyenv())
gc_get_lib <- function() {
  if (is.null(gc_lib$fns)) {
    env <- new.env(parent = globalenv())
    tryCatch(Rcpp::sourceCpp(code = gc_index_cpp, env = env), error = function(err) {
      stop("GC matching needs a C++ compiler for Rcpp, and compiling failed:\n", conditionMessage(err), call. = FALSE)
    })
    gc_lib$fns <- env
  }
  gc_lib$fns
}
# window indexes, shared by all randomizers in one R session
gc_index_store <- new.env(parent = emptyenv())

# set up GC matching
prepare_gc <- function(genome, tolerance = 0.02, flank = 0L, max_n = 0.1) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0) {
    stop("GC tolerance must be a single number > 0, got: ", paste(tolerance, collapse = ", "))
  }
  if (!is.numeric(flank) || length(flank) != 1L || flank < 0) {
    stop("GC flank must be a single number >= 0, got: ", paste(flank, collapse = ", "))
  }
  if (!is.numeric(max_n) || length(max_n) != 1L || max_n < 0 || max_n > 1) {
    stop("GC max_n must be between 0 and 1, got: ", paste(max_n, collapse = ", "))
  }

  # needs a BSgenome
  sequences <- if (inherits(genome, "BSgenome") || methods::is(genome, "DNAStringSet")) {
    genome
  } else {
    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("GC matching requires the BSgenome package.")
    }
    tryCatch(
      BSgenome::getBSgenome(genome, masked = FALSE),
      error = function(err) {
        stop(
          "GC matching could not load a BSgenome for '", genome, "'.\n",
          conditionMessage(err), "\n",
          "Run: nullsens genomes install ", genome
        )
      }
    )
  }

  genome_seqs <- if (inherits(sequences, "BSgenome")) {
    GenomeInfoDb::seqlengths(sequences)
  } else {
    stats::setNames(as.numeric(Biostrings::width(sequences)), names(sequences))
  }

  # index for one chromosome, built on first use
  index_cache <- new.env(parent = emptyenv())
  source_name <- if (is.character(genome)) genome else NULL
  fetch_index <- function(chr) {
    key <- if (is.null(source_name)) NULL else paste(source_name, chr, sep = "|")
    hit <- if (is.null(key)) index_cache[[chr]] else gc_index_store[[key]]
    if (is.null(hit)) {
      seq <- sequences[[chr]]
      if (methods::is(seq, "MaskedDNAString")) seq <- Biostrings::unmasked(seq)
      hit <- gc_get_lib()$gc_build_index(charToRaw(as.character(seq)))
      if (is.null(key)) index_cache[[chr]] <- hit else gc_index_store[[key]] <- hit
    }
    hit
  }

  list(
    fetch_index = fetch_index,
    orig_cache  = new.env(parent = emptyenv()),   # cached GC of the originals
    seqnames    = names(genome_seqs),
    seq_lengths = genome_seqs,
    source      = if (is.character(genome)) genome else class(genome)[1],
    tolerance   = tolerance,
    flank       = as.integer(flank),
    max_n       = max_n,
    n_sequences = length(genome_seqs)
  )
}

# GC and N fraction of each region (with flank)
gc_content <- function(context, regions) {
  n <- length(regions)
  gc <- rep(NA_real_, n)
  n_frac <- rep(1, n)
  if (n == 0L) return(list(gc = gc, n_frac = n_frac))

  chroms <- as.character(GenomicRanges::seqnames(regions))
  chrom_len <- unname(context$seq_lengths[chroms])
  left <- pmax(1, GenomicRanges::start(regions) - context$flank)
  right <- pmin(chrom_len, GenomicRanges::end(regions) + context$flank)
  # skip windows that are off the chromosome
  usable <- !is.na(chrom_len) & left <= right
  if (!any(usable)) return(list(gc = gc, n_frac = n_frac))

  uchr <- unique(chroms[usable])
  counts <- gc_get_lib()$gc_window_stats(lapply(uchr, context$fetch_index), match(chroms[usable], uchr),
                                         as.integer(left[usable]), as.integer(right[usable]))
  widths <- right[usable] - left[usable] + 1L
  called <- widths - counts$n
  gc[usable] <- ifelse(called > 0, counts$gc / called, NA_real_)
  n_frac[usable] <- counts$n / widths
  list(gc = gc, n_frac = n_frac)
}

gc_region_keys <- function(regions) {
  paste(
    as.character(GenomicRanges::seqnames(regions)),
    GenomicRanges::start(regions), GenomicRanges::end(regions), sep = ":"
  )
}

# like gc_content(), but each original region is computed only once
gc_content_original <- function(context, regions) {
  cache <- context$orig_cache
  if (is.null(cache) || length(regions) == 0L) {
    return(gc_content(context, regions))
  }
  keys <- gc_region_keys(regions)
  n_keys <- paste0("n|", keys)
  hit <- mget(keys, envir = cache, ifnotfound = list(NULL))
  miss <- vapply(hit, is.null, logical(1))
  gc <- numeric(length(keys))
  n_frac <- numeric(length(keys))
  if (any(!miss)) {
    gc[!miss] <- unlist(hit[!miss], use.names = FALSE)
    n_frac[!miss] <- unlist(mget(n_keys[!miss], envir = cache), use.names = FALSE)
  }
  if (any(miss)) {
    fresh <- gc_content(context, regions[miss])
    gc[miss] <- fresh$gc
    n_frac[miss] <- fresh$n_frac
    first <- !duplicated(keys[miss])
    list2env(stats::setNames(as.list(fresh$gc[first]), keys[miss][first]), envir = cache)
    list2env(stats::setNames(as.list(fresh$n_frac[first]), n_keys[miss][first]), envir = cache)
  }
  list(gc = gc, n_frac = n_frac)
}

# original_gc: GC of each original, if the caller already has it (build_randomizer does)
validate_gc <- function(context, original, candidate, original_gc = NULL) {
  if (length(original) != length(candidate)) {
    stop("validate_gc: original and candidate must have the same length and order.")
  }
  if (!is.null(original_gc) && length(original_gc) != length(original)) {
    stop("validate_gc: original_gc must have one value per original region.")
  }

  orig <- list(gc = if (is.null(original_gc)) gc_content_original(context, original)$gc else original_gc)
  cand <- gc_content(context, candidate)

  valid <- abs(cand$gc - orig$gc) <= context$tolerance & cand$n_frac <= context$max_n

  candidate_chroms <- as.character(GenomicRanges::seqnames(candidate))
  candidate_lengths <- unname(context$seq_lengths[candidate_chroms])
  in_bounds <- !is.na(candidate_lengths) &
    GenomicRanges::start(candidate) >= 1L &
    GenomicRanges::end(candidate) <= candidate_lengths
  valid[!in_bounds |
          GenomicRanges::width(candidate) != GenomicRanges::width(original)] <- FALSE

  # original with no GC: only the N check applies
  no_gc <- is.na(orig$gc)
  valid[no_gc] <- cand$n_frac[no_gc] <= context$max_n

  valid[is.na(valid)] <- FALSE
  valid
}

# ---- Direct samplers (not used by the current build_randomizer) ----
# Draw each region straight from windows that pass validate_gc, instead of propose-and-retry.
#   make_gc_exact_sampler(context, proposal_genome, mask, cores)  -> function(region): one region, or NULL
#   make_gc_proposal_sampler(proposal_genome, mask)               -> function(A): one window per region
# Build the exact sampler first.
gc_registry <- new.env(parent = emptyenv())
gc_registry$entries <- list()

gc_find_entry <- function(proposal_genome, mask) {
  for (candidate in rev(gc_registry$entries)) {
    if (identical(candidate$genome, proposal_genome) && identical(candidate$mask, mask)) {
      return(candidate)
    }
  }
  NULL
}

# one GC-compatible window per region (a region with none is returned unchanged)
make_gc_proposal_sampler <- function(proposal_genome, mask = NULL) {
  function(A) {
    if (length(A) == 0L) return(A)
    entry <- gc_find_entry(proposal_genome, mask)
    if (is.null(entry)) {
      stop("make_gc_exact_sampler() must be called with the same genome and mask before the proposal sampler is used.")
    }
    drawn <- entry$draw_many(A)
    GenomicRanges::GRanges(
      seqnames = ifelse(drawn$found, drawn$chromosome, as.character(GenomicRanges::seqnames(A))),
      ranges = IRanges::IRanges(
        start = ifelse(drawn$found, drawn$start, GenomicRanges::start(A)),
        width = GenomicRanges::width(A)
      ),
      seqinfo = GenomeInfoDb::seqinfo(proposal_genome)
    )
  }
}

# cores is only checked, not used
make_gc_exact_sampler <- function(context, proposal_genome, mask = NULL, cores = 1L) {
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1L || cores != as.integer(cores)) {
    stop("cores must be a positive integer.")
  }
  lib <- gc_get_lib()

  chroms <- unique(as.character(GenomicRanges::seqnames(proposal_genome)))
  chroms <- chroms[chroms %in% context$seqnames]
  genome_ranges <- stats::setNames(lapply(chroms, function(chrom) {
    keep <- as.character(GenomicRanges::seqnames(proposal_genome)) == chrom
    IRanges::reduce(GenomicRanges::ranges(proposal_genome[keep]))
  }), chroms)
  mask_ranges <- stats::setNames(lapply(chroms, function(chrom) {
    if (is.null(mask)) return(IRanges::IRanges())
    keep <- as.character(GenomicRanges::seqnames(mask)) == chrom
    IRanges::reduce(GenomicRanges::ranges(mask[keep]))
  }), chroms)
  has_mask <- any(lengths(mask_ranges) > 0L)

  indexes <- NULL
  load_indexes <- function() {
    if (is.null(indexes)) indexes <<- lapply(chroms, context$fetch_index)
  }

  # all legal starts for each width, with cumulative counts for uniform draws
  legal_cache <- new.env(parent = emptyenv())
  get_legal_index <- function(w) {
    key <- as.character(w)
    ix <- legal_cache[[key]]
    if (!is.null(ix)) return(ix)
    st <- list(); ch <- list(); wd <- list()
    for (ci in seq_along(chroms)) {
      gr <- genome_ranges[[ci]]
      first <- IRanges::start(gr)
      last <- IRanges::end(gr) - w + 1L
      keep <- last >= first
      if (any(keep)) {
        st[[length(st) + 1L]] <- first[keep]
        wd[[length(wd) + 1L]] <- last[keep] - first[keep] + 1L
        ch[[length(ch) + 1L]] <- rep(ci, sum(keep))
      }
    }
    ix <- if (length(st) == 0L) {
      list(total = 0)
    } else {
      wd <- unlist(wd)
      list(starts = unlist(st), chrom = unlist(ch), cumulative = cumsum(as.double(wd)), total = sum(as.double(wd)))
    }
    legal_cache[[key]] <- ix
    ix
  }

  drop_masked <- function(chrom_i, starts, width) {
    keep <- rep(TRUE, length(starts))
    if (has_mask && length(starts) > 0L) {
      for (ci in unique(chrom_i)) {
        masks <- mask_ranges[[ci]]
        if (length(masks) == 0L) next
        here <- which(chrom_i == ci)
        keep[here] <- !IRanges::overlapsAny(
          IRanges::IRanges(start = starts[here], width = if (length(width) == 1L) width else width[here]), masks)
      }
    }
    keep
  }

  # try random legal starts (16 per region, then 8x more each round); keep the first that passes the GC test and avoids the mask
  fast_draw <- function(widths, target) {
    n <- length(widths)
    chrom_i <- integer(n); start <- integer(n); found <- logical(n)
    pending <- seq_len(n)
    tries <- 16L
    while (length(pending) > 0L && tries <= 524288L) {
      per_chunk <- max(1L, as.integer(floor(3e6 / tries)))
      for (chunk in split(pending, (seq_along(pending) - 1L) %/% per_chunk)) {
        reg <- list(); pos_start <- list(); pos_chrom <- list()
        for (group in split(chunk, widths[chunk])) {
          ix <- get_legal_index(widths[group[1L]])
          if (ix$total == 0) next
          r <- rep(group, each = tries)
          draws <- floor(stats::runif(length(r)) * ix$total)
          ri <- findInterval(draws, ix$cumulative) + 1L
          prev <- numeric(length(ri))
          has_prev <- ri > 1L
          prev[has_prev] <- ix$cumulative[ri[has_prev] - 1L]
          reg[[length(reg) + 1L]] <- r
          pos_start[[length(pos_start) + 1L]] <- ix$starts[ri] + as.integer(draws - prev)
          pos_chrom[[length(pos_chrom) + 1L]] <- ix$chrom[ri]
        }
        if (length(reg) == 0L) next
        reg <- unlist(reg, use.names = FALSE)
        pos_start <- unlist(pos_start, use.names = FALSE)
        pos_chrom <- unlist(pos_chrom, use.names = FALSE)
        w <- widths[reg]
        ok <- lib$gc_window_ok(indexes, pos_chrom, pos_start, w, target[reg],
                               context$tolerance, context$max_n, context$flank)
        survivors <- which(ok)
        if (length(survivors) > 0L) {
          ok[survivors] <- drop_masked(pos_chrom[survivors], pos_start[survivors], w[survivors])
        }
        hit <- which(ok)
        if (length(hit) > 0L) {
          first_hit <- hit[!duplicated(reg[hit])]
          r <- reg[first_hit]
          found[r] <- TRUE
          start[r] <- pos_start[first_hit]
          chrom_i[r] <- pos_chrom[first_hit]
        }
      }
      pending <- pending[!found[pending]]
      tries <- tries * 8L
    }
    list(chrom = chrom_i, start = start, found = found)
  }

  # rare targets: list every compatible window once and draw from that list
  rare <- new.env(parent = emptyenv())
  rare_keys <- character(0)
  target_key <- function(width, target) {
    paste(width, ifelse(is.na(target), "NA", sprintf("%.17g", target)), sep = ":")
  }
  rare_pool <- function(width, target) {
    key <- target_key(width, target)
    pool <- rare[[key]]
    if (is.null(pool)) {
      scan_chromosome <- function(ci) {
        gr <- genome_ranges[[ci]]
        starts <- integer(0)
        for (k in seq_along(gr)) {
          first <- IRanges::start(gr)[k]
          last <- IRanges::end(gr)[k] - width + 1L
          if (last < first) next
          found <- lib$gc_all_compatible_starts(indexes, ci, first, last, width, target,
                                                context$tolerance, context$max_n, context$flank, 1000000L)
          if (found$overflow) stop("too many compatible windows to list for width ", width, " bp")
          starts <- c(starts, found$starts)
        }
        starts
      }
      per_chrom <- lapply(seq_along(chroms), scan_chromosome)
      chrom_i <- rep(seq_along(chroms), lengths(per_chrom))
      starts <- unlist(per_chrom, use.names = FALSE)
      keep <- drop_masked(chrom_i, starts, width)
      pool <- list(chrom = chrom_i[keep], start = starts[keep])
      rare[[key]] <- pool
    }
    pool
  }

  # one uniform draw per region
  draw_many <- function(original) {
    n <- length(original)
    chromosome <- character(n); start <- integer(n); found <- logical(n)
    if (n == 0L) return(list(chromosome = chromosome, start = start, found = found))
    load_indexes()
    widths <- as.integer(GenomicRanges::width(original))
    target <- gc_content_original(context, original)$gc
    key <- target_key(widths, target)

    try_fast <- which(!(key %in% rare_keys))
    if (length(try_fast) > 0L) {
      res <- fast_draw(widths[try_fast], target[try_fast])
      f <- which(res$found)
      found[try_fast[f]] <- TRUE
      start[try_fast[f]] <- res$start[f]
      chromosome[try_fast[f]] <- chroms[res$chrom[f]]
    }
    for (i in which(!found)) {
      rare_keys <<- union(rare_keys, key[i])
      pool <- rare_pool(widths[i], target[i])
      m <- length(pool$start)
      if (m == 0L) next
      pick <- if (m == 1L) 1L else sample.int(m, 1L)
      chromosome[i] <- chroms[pool$chrom[pick]]
      start[i] <- pool$start[pick]
      found[i] <- TRUE
    }
    list(chromosome = chromosome, start = start, found = found)
  }

  # one region at a time; NULL if there is no compatible window
  sampler <- function(original) {
    if (length(original) != 1L) {
      stop("The exact GC sampler accepts one original region at a time.")
    }
    drawn <- draw_many(original)
    if (!drawn$found[1L]) {
      return(NULL)
    }
    GenomicRanges::GRanges(
      seqnames = drawn$chromosome,
      ranges = IRanges::IRanges(start = drawn$start, width = GenomicRanges::width(original)),
      seqinfo = GenomeInfoDb::seqinfo(proposal_genome)
    )
  }

  gc_registry$entries <- c(
    utils::tail(gc_registry$entries, 7L),
    list(list(genome = proposal_genome, mask = mask, draw_many = draw_many))
  )
  sampler
}
