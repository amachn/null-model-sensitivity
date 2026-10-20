# ~~~ src/nullsens/r/regioner_runner.R ~~~
# - R backend entry point
#   - takes JSON input from Python
#   - validation?
#   - loads the region files (rtracklayer)
#   - builds the randomizer (outsource to another fn/file for abstraction?)
#   - runs regioneR::permTest once
#   - outputs the result as JSON
suppressPackageStartupMessages({
  library(jsonlite)
})

source("nulls/length.R")

# read the JSON request from stdin
input <- paste(readLines(file("stdin"), warn = FALSE), collapse = "\n")
config <- fromJSON(input)
