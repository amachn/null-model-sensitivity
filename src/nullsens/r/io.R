# ~~~ src/nullsens/r/io.R ~~~
# - imports genomic regions into a common GRanges representation

library(rtracklayer)

get_file_extension <- function(path) {
  base <- tolower(basename(path))
  base <- sub("\\.(gz|bgz|bz2|xz|zip)$", "", base)
  tools::file_ext(base)
}

find_column <- function(data, candidates) {
  cols <- names(data)
  hits <- match(tolower(candidates), tolower(cols))
  hits <- hits[!is.na(hits)]
  
  if (length(hits) == 0) {
    stop(
      "None of the expected columns found: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(cols, collapse = ", ")
    )
  }
  
  cols[hits[1]]
}

read_tabular <- function(path) {
  extension <- get_file_extension(path)
  
  if (extension %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel files: ", path)
    }
    return(as.data.frame(readxl::read_excel(path)))
  }
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    sep <- if (extension == "csv") "," else "\t"
    return(utils::read.table(path, header = TRUE, sep = sep,
                             stringsAsFactors = FALSE, check.names = FALSE))
  }
  as.data.frame(data.table::fread(path, data.table = FALSE))
}

import_tabular_regions <- function(path, zero_based = FALSE) {
  
  data <- read_tabular(path)
  
  start_col <- find_column(
    data, c("start", "st", "begin", "chromstart")
  )
  
  end_col <- find_column(
    data, c("end", "chromend")
  )
  
  chr_col <- find_column(
    data, c("chr", "chromosome", "chrom", "seqnames")
  )
  
  gr <- makeGRangesFromDataFrame(
    data,
    keep.extra.columns         = FALSE,
    seqnames.field             = chr_col,
    start.field                = start_col,
    end.field                  = end_col,
    starts.in.df.are.0based    = zero_based
  )
  
  gr
}

import_bigbed_regions <- function(path) {
  
  bb_file <- BigBedFile(path)
  gr <- import(bb_file)
  
  mcols(gr) <- NULL
  gr
}


import_regions <- function(path) {
  
  if (!file.exists(path)) {
    stop("Region file not found: ", path)
  }
  
  extension <- get_file_extension(path)
  
  if (extension %in% c("txt", "tsv", "csv", "xlsx", "xls")) {
    return(import_tabular_regions(path))
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
