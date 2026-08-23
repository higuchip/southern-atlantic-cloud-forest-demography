required_inputs <- c(
  "data/tree_census.csv",
  "data/environmental.csv"
)

if (!all(file.exists(required_inputs))) {
  stop("Run this script from the repository root; one or more inputs are missing.")
}

source("R/01_run_analysis.R", encoding = "UTF-8")
source("R/02_export_paper_table.R", encoding = "UTF-8")

if (file.exists("R/03_build_paper_figures.R")) {
  source("R/03_build_paper_figures.R", encoding = "UTF-8")
}

cat("Analysis, table, and figures 2-3 finished successfully.\n")
