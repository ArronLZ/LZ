#' delete the duplicate name of the first column of dataframe
#'
#' @param eset dataframe, the first column must be gene, type: character
#' @param col.by character, delete the duplicate by this column
#' @param col.del character, manually delete some column
#' @param auto.del.character logical, if auto delete the column those class is character
#'
#' @return dataframe
#' @export
#' @import dplyr
#' @importFrom tibble column_to_rownames
#' 
#' @author Jiang
#' @examples
#' \dontrun{
#' quchong(eset, col.by="id")
#' head(eset)[,1:3]
#'           gene_id cont-1 cont-2
#'   ENSG00000186827      0      1
#'   ENSG00000186891     63     50
#'   ENSG00000160072   1218   1023
#' }
quchong <- function(eset, col.by, col.del=NULL, auto.del.character=T) {
  if (!is.null(col.del)) {
    eset <- eset[, -col.del]
    cat("!!! 请注意 ", col.del, "列被删除\n")
  }
  if (!col.by %in% names(eset)) {
    stop("col.by必须为eset的列名之一，请重新设置col.by参数")
  }
  if (sum(sapply(eset, is.numeric)) < 1) {
    stop("eset数据整体格式不对，数据中没有一列是数值型变量，请正确载入数据")
  }
  # col.by.num <- match(col.by, names(eset))
  if (is.numeric(eset[, col.by]) == T) {
    stop("col.by列的变量类型必须为character, 请注意数据col.by列是否为gene")
  }
  # names(eset)[col.by.num] <- "IDtemplz"
  eset$MEAN <- abs(rowMeans(eset[, sapply(eset, is.numeric)]))
  eset <- eset[!duplicated(eset[, col.by]), ]
  rownames(eset) <- NULL
  eset <- tibble::column_to_rownames(eset, var = col.by)
  if (!is.null(eset$MEAN)) { eset$MEAN <- NULL }
  if (auto.del.character & (sum(!sapply(eset, is.numeric)) > 0) ) {
    eset <- eset[, !sapply(eset, is.numeric)]
    cat("!!! 请注意 ", names(eset)[!sapply(eset, is.numeric)], "列被删除\n")
  }
  return(eset)
}


#' clean a dir name
#' @description clean a dir name
#' @param dir character the dir name
#'
#' @return dirname
#' @export
#' 
#' @author Jiang
dirclean <- function(dir) {
  dir <- trimws(dir, which = "both", whitespace = "[ \t\r\n\\/\\\\]")
  return(dir)
}


#' create dir
#' @description mkdir, in deg analysis usually used to make output dir
#'
#' @param dir character, the dir name
#'
#' @return # if the dir no exist the function will create the dir
#' @export
#'
#' @author Jiang
mkdir <- function(dir) {
  if(dir.exists(dir)) {
    cat(paste0("Warning: The directory ", dir, "already exists, do you want to 
    delete the entire directory and re-create it ? ! 
    IF YOU ARE SURE WHAT YOU ARE DOING, please use the function mkdir.p.
    Otherwise, you can ignore this warning message.\n"))
  } else {
    dir.create(dir, recursive = T, showWarnings = T) 
  }
}


#' Create a folder, regardless of whether the folder already exists
#' @description Create a folder, regardless of whether the folder already exists, overwrite if it exists, in deg analysis usually used to make output dir.
#'
#' @param dir character, the dir name
#'
#' @return # Create a folder regardless of whether the folder already exists.
#' @export
#'
#' @examples #
mkdir.p <- function(dir) {
  unlink(dir, recursive = T)
  dir.create(dir, recursive = T, showWarnings = T)
}

#' find the NA(NA, grepl("\\s|^$", lie))
#' @description find the NA(NA, grepl("\\s|^$", lie))
#' @param lie vector the column or vector of your data
#'
#' @return vector
#' @export
#'
#' @examples # lz_isna()
lz_isna <- function(lie) {
  # 查空值
  re <- is.na(lie) | grepl("\\s|^$", lie)
  return(re)
}


#' show part info of a dataframe
#' @description show part info of a dataframe
#' @param df dataframe data
#' @param ncol number default is 4, the end colomn index
#' @param nrow number default is 3, the end colomn index
#' @param rowstart number default is 1, the start row index
#' @param colstart number default is 1, the start colomn index
#' 
#' @return #
#' @export
#'
#' @examples # df %>% lzhead
lzhead <- function(df, ncol = 4, colstart = 1) {
  if (ncol(df) > 4) { 
    head(df)[, colstart:ncol] 
  } else {
    head(df)
  }
}


#' fpkm to tmp
#' @description fpkm to tmp
#' @param fpkm numeric the data should be no log data
#'
#' @return tpm
#' @export
#'
#' @examples # fpkmToTpm()
fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


#' Annot Data -- gencode.v22.annot
#' @description Annot Data -- gencode.v22.annot
#' @param ... character, the alias of bulit-in dataset, ex. gencode.v22.annot
#' keggpathway.gmt
#' 
#' @return list Large list
#' @export
#' 
#' @author Jiang
#' @examples
#' #  data(gencode.v22.annot, package = "LZ")
#' #  data(keggpathway.gmt, package = "LZ")
getData <- function(..., package = "LZ") {
  data(..., package = package)
}
