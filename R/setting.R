#' install package by using CRAN or bioconductor
#' @description install package by using CRAN or bioconductor, if the package is not in
#' the CRAN, then the function will use bioconductor method to install the package.
#'
#' @param package charater the package name you want to install
#'
#' @return #
#' @export
#' @importFrom BiocManager install
#'
#' @examples #
#' # packs <- c('dplyr', 'stringr', 'DESeq2')
#' # for (p in packs) {
#' #   LZ::install(p)
#' # }
#' # or
#' # LZ::install("ggrepel")
#' @author Jiang
install <- function(package) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
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


#' check if the package is been installed
#' @description check if the package is been installed
#'
#' @param packname character the R package name
#'
#' @return # message in console
#' @export
#'
#' @author Jiang
require_pack <- function(packname) {
  if (!requireNamespace(packname, quietly = TRUE)) {
    cat(paste("Error: '", packname, "' package is required but not installed.\n", sep = ""))
    cat(paste("Please install the package using install.packages('", packname, "')\n",sep = ""))
  }
}


#' set http proxy for R session
#' @description set http proxy for R session
#' @param ngate the proxy address, for example, if you use clash, the address is default "http://127.0.0.1:7890"
#'
#' @return no return，use unsetproxy to cancel the proxy
#' @export
#'
#' @examples
#' # setproxy()
#' # setproxy("http://127.0.0.1:7890")
#' @author Jiang
setproxy <- function(ngate="http://127.0.0.1:7890") {
  Sys.setenv('http_proxy'=ngate)
  Sys.setenv('https_proxy'=ngate)
}


#' install package by using CRAN or bioconductor
#' @description install package by using CRAN or bioconductor, if the package is not in
#' the CRAN, then the function will use bioconductor method to install the package.
#' @param packs charater(n), the package name you want to install
#'
#' @return no
#' @export
#'
#' @author Jiang
p_install <- function(packs) {
  for (p in packs) {
    install(p)
  }
}
