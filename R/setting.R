#' install package by using CRAN or bioconductor
#' @description install package by using CRAN or bioconductor, if the package is not in
#' the CRAN, then the function will use bioconductor method to install the package.
#'
#' @param package charater the package name you want to install
#'
#' @return #
#' @export
#'
#' @examples #
#' # packs <- c('dplyr', 'stringr', 'DESeq2')
#' # for (p in packs) {
#' #   LZ::install(p)
#' # }
#' # or
#' # LZ::install("ggrepel")
install <- function(package) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install(version = "3.17", ask = F)
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    cat("Attempting to install", package, "using install.packages()\n")
    install.packages(package)
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    cat("Installation using install.packages() failed. Attempting to install",
        package, "using BiocManager::install()\n")
    BiocManager::install(package, update = F, ask = F)
  }
  if (requireNamespace(package, quietly = TRUE)) {
    cat("Install Complete:  ", package, "\n")
  }
}
