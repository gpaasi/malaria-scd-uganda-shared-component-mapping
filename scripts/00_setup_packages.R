suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(exactextractr)
  library(gdistance)
  library(future.apply)
})

# global options
options(stringsAsFactors = FALSE)

# helper to check required packages are loaded
required <- c("data.table","sf","terra","exactextractr","gdistance","future.apply")
loaded <- sapply(required, function(pkg) requireNamespace(pkg, quietly=TRUE))
if (!all(loaded)) {
  missing <- required[!loaded]
  warning("The following packages are not installed: ", paste(missing, collapse = ", "), 
          ". Install with install.packages(", paste0('"', missing, '"', collapse = ", "), ").")
}
