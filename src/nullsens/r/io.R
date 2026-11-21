# ~~~ src/nullsens/r/io.R ~~~
# - imports genomic region inputs into a common GRanges representation
# - normalizes imported regions to the requested reference genome

# normalize column names for flexible matching
normalize_column_name <- function(name) {
  tolower(gsub("[^a-z0-9]", "", name))
}

# get the file extension, ignoring compression extensions
get_file_extension <- function(path) {
  base <- tolower(basename(path))
  base <- sub("\\.(gz|bgz|bz2|xz|zip)$", "", base)
  tools::file_ext(base)
}

# find the first matching column from a set of common names
find_column <- function(data, candidates) {
  columns <- names(data)
  
  normalized_columns <- normalize_column_name(columns)
  normalized_candidates <- normalize_column_name(candidates)
  
  match_index <- match(normalized_columns, normalized_candidates, nomatch = 0L)
  match_index <- match_index[match_index > 0L]
  
  if (length(match_index) == 0L) { return(NULL) }
  columns[match_index[1]]
}

# read a generic tabular file into a data.frame
read_tabular <- function(path, header = TRUE) {
  extension <- get_file_extension(path)
  
  if (extension %in% c("xlsx", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel files: ", path)
    }
    return(as.data.frame(readxl::read_excel(path, col_names = header)))
  }
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    sep <- if (extension == "csv") "," else "\t"
    return(utils::read.table(path, header = header, sep = sep,
                             stringsAsFactors = FALSE, check.names = FALSE))
  }
  
  as.data.frame(data.table::fread(path, header = header, data.table = FALSE))
}

# import a generic tabular genomic-region file
import_tabular_regions <- function(path, zero_based = NULL) {
  data <- read_tabular(path, header = TRUE)
  
  chr_col <- find_column(data, c("chr", "chrom", "chromosome", "seqnames"))
  start_col <- find_column(data, c("st", "start", "chromStart", "begin"))
  end_col <- find_column(data, c("end", "stop", "chromEnd"))
  
  # if coord cols aren't found, retry as a headerless BED-like file
  if (is.null(chr_col) || is.null(start_col) || is.null(end_col)) {
    data <- read_tabular(path, header = FALSE)
    if (ncol(data) < 3L) {
      stop("Unable to identify genomic coordinates in: ", path)
    }
    
    names(data)[1:3] <- c("chrom", "start", "end")
    chr_col <- "chrom"
    start_col <- "start"
    end_col <- "end"
  }
  
  # infer coordinate convention when not specified
  # * chromStart is the standard BED field name, and thus 0-based
  if (is.null(zero_based)) {
    zero_based <- normalize_column_name(start_col) == "chromstart"
  }
  
  GenomicRanges::makeGRangesFromDataFrame(
    data,
    seqnames.field          = chr_col,
    start.field             = start_col,
    end.field               = end_col,
    starts.in.df.are.0based = zero_based,
    keep.extra.columns      = FALSE
  )
}

# import a BigBed region file
# * rtracklayer's BigBed impl is not supported on Windows
# TODO: add BigBed import support on Windows using UCSC bigBedToBed?
import_bigbed_regions <- function(path) {
  bb_file <- rtracklayer::BigBedFile(path)
  gr <- rtracklayer::import(bb_file)
  S4Vectors::mcols(gr) <- NULL
  gr
}

# import an analysis region set into GRanges
import_regions <- function(path, zero_based = NULL) {
  if (!file.exists(path)) {
    stop("Region file not found: ", path)
  }
  
  extension <- get_file_extension(path)
  
  if (extension %in% c("txt", "tsv", "csv", "xlsx", "xls")) {
    return(import_tabular_regions(path, zero_based))
  }
  
  if (extension %in% c("bigbed", "bb")) {
    return(import_bigbed_regions(path))
  }
  
  tryCatch(
    rtracklayer::import(path),
    error = function(err) {
      stop("Unable to import region file: ", path, "\n", conditionMessage(err))
    }
  )
}

# normalize a genome spec into a GRanges object
# caller provides either:
# - genome name understood by regioneR, e.g. "hg38"
# - an existing genome GRanges object
as_genome_gr <- function(genome) {
  genome_gr <- if (methods::is(genome, "GRanges")) {
    genome 
  } else {
    suppressWarnings(
      suppressMessages(
        regioneR::getGenome(genome)
      )
    )
  }
  
  GenomicRanges::strand(genome_gr) <- "*"
  
  if (anyNA(GenomeInfoDb::seqlengths(genome_gr))) {
    # fall back to the range ends as chromosome lengths
    GenomeInfoDb::seqlengths(genome_gr) <- as.integer(
      tapply(GenomicRanges::end(genome_gr), 
             as.character(GenomicRanges::seqnames(genome_gr)),
             max)[GenomeInfoDb::seqlevels(genome_gr)]
    )
  }
  
  genome_gr
}

# conform imported genomic regions to the reference genome
conform_regions <- function(regions, genome_gr) {
  region_seqs <- GenomeInfoDb::seqlevels(regions)
  genome_seqs <- GenomeInfoDb::seqlevels(genome_gr)
  
  # if no direct chr-name overlap, try changing the region set to the
  # naming style used by the selected genome
  if (!any(region_seqs %in% genome_seqs)) {
    target_style <- tryCatch(
      GenomeInfoDb::seqlevelsStyle(genome_gr)[1],
      error = function(err) NA_character_
    )
    
    if (!is.na(target_style)) {
      tryCatch(
        GenomeInfoDb::seqlevelsStyle(regions) <- target_style,
        error = function(err) NULL
      )
    }
    
    region_seqs <- GenomeInfoDb::seqlevels(regions)
  }
  
  common <- intersect(region_seqs, genome_seqs)
  if (length(common) == 0L) {
    stop("No chromosome names in the region set match the requested genome.")
  }
  
  original_n <- length(regions)
  regions <- GenomeInfoDb::keepSeqlevels(
    regions, common, pruning.mode = "coarse")
  
  if (length(regions) < original_n) {
    warning(
      original_n - length(regions), 
      "region(s) were droped because their chromosomes aren't present in",
      " the requested genome."
    )
  }
  
  # assign the chr lengths + other Seqinfo metadata from the selected genome
  GenomeInfoDb::seqinfo(regions) <-
    GenomeInfoDb::seqinfo(genome_gr)[
      GenomeInfoDb::seqlevels(regions)
    ]
  
  # clamp any ranges outside known chr boundaries
  regions <- suppressWarnings(GenomicRanges::trim(regions))
  
  GenomicRanges::strand(regions) <- "*"
  
  regions
}
