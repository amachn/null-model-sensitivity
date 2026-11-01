# ~~~ src/nullsens/r/io.R ~~~
# - imports genomic regions into a common GRanges representation

get_file_extension <- function(path) {}
find_column <- function(data, candidates) {
  
}

import_tabular_regions <- function(path) {
  start_col <- find_column(
    data, c("start", "st", "begin", "chromstart")
  )
}
import_bigbed_regions <- function(path) {}


# import an analysis region set into GRanges
import_regions <- function(path) {
  extension <- get_file_extension(path)
  
  if (extension %in% c("txt", "tsv")) {
    return(import_tabular_regions(path))
  }
  
  if (extension %in% c("bigbed")) {
    return(import_bigbed_regions(path))
  }
  
  tryCatch(
    rtracklayer::import(path),
    error = function(err) {
      stop("Unable to import region file: ", path, "\n", conditionMessage(err))
    }
  )
}