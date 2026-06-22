#!/usr/bin/env Rscript

# Dependency installer/checker for NanoTel.R
# Usage: Rscript requirements.R

minimum_r_version <- "4.2.2"

if (getRversion() < minimum_r_version) {
  stop(
    sprintf(
      "NanoTel.R requires R >= %s (installed: %s).",
      minimum_r_version,
      getRversion()
    ),
    call. = FALSE
  )
}

cran_packages <- c(
  "optparse",
  "conflicted",
  "tidyverse",
  "logr",
  "future",
  "ggprism",
  "testit",
  "survival"
)

bioconductor_packages <- c(
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "Biostrings"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0L) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0L) {
  message(
    "Installing missing Bioconductor packages: ",
    paste(missing_bioc, collapse = ", ")
  )
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

all_packages <- c(cran_packages, bioconductor_packages)
still_missing <- all_packages[
  !vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0L) {
  stop(
    "The following dependencies could not be installed: ",
    paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}

message("All NanoTel.R package requirements are installed.")
message("R version: ", getRversion())
message("Operating system: ", Sys.info()[["sysname"]])

if (.Platform$OS.type == "windows") {
  message(
    paste(
      "Windows note: future::multicore, FORK clusters, and multicore mclapply",
      "are not supported. NanoTel.R must use a multisession/PSOCK backend."
    )
  )
} else {
  message(
    paste(
      "Unix note: the script can use future::multicore and FORK workers,",
      "subject to the safety restrictions of the R frontend being used."
    )
  )
}
