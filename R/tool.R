#' set http proxy for R session
#' @description set http proxy for R session
#'
#' @param ngate the proxy address, for example, if you use clash, the address is default "http://127.0.0.1:7890"
#'
#' @return no return，use unsetproxy to cancel the proxy
#' @export
#'
#' @examples #
setproxy <- function(ngate="http://127.0.0.1:7890") {
  Sys.setenv('http_proxy'=ngate)
  Sys.setenv('https_proxy'=ngate)
}


#' delete the duplicate name of the select coloumn of dataframe
#' delete the duplicate name of the select coloumn of dataframe
#' @param eset dataframe, the first coloumn must be gene, type: character
#'
#' @return no return
#' @export
#' @importFrom tibble column_to_rownames
#'
#' @examples #
#' #           gene_id cont-1 cont-2
#' #   ENSG00000186827      0      1
#' #   ENSG00000186891     63     50
#' #   ENSG00000160072   1218   1023
quchong <- function(eset) {
  if (is.numeric(eset[,1]) == T) {
    stop("第一列必须为character, 请注意数据第一列是否为gene")
  }
  names(eset)[1] <- "ID"
  eset$MEAN <- abs(rowMeans(eset[,2:ncol(eset)]))
  eset <- eset %>% arrange(ID, desc(MEAN))
  eset <- eset[!duplicated(eset$ID), ]
  rownames(eset) <- NULL
  eset <- tibble::column_to_rownames(eset, var = "ID")
  if (!is.null(eset$MEAN)) { eset$MEAN <- NULL }
  return(eset)
}


#' check if the package is been installed
#' @description check if the package is been installed
#'
#' @param packname character the R package name
#'
#' @return # message in console
#' @export
#'
#' @examples #
require_pack <- function(packname) {
  if (!requireNamespace(packname, quietly = TRUE)) {
    cat(paste("Error: '", packname, "' package is required but not installed.\n", sep = ""))
    cat(paste("Please install the package using install.packages('", packname, "')\n",sep = ""))
  }
}


#' create dir
#' @description mkdir, in deg analysis usually used to make output dir.
#'
#' @param dir character, the dir name
#'
#' @return # if the dir no exist the function will create the dir.
#' @export
#'
#' @examples #
mkdir <- function(dir) {
  if(!dir.exists(dir)) { dir.create(dir, recursive = T, showWarnings = T) }
}
