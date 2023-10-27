#' delete the duplicate name of the first column of dataframe
#' @param eset dataframe, the first column must be gene, type: character
#'
#' @return dataframe
#' @export
#' @import dplyr
#' @importFrom tibble column_to_rownames
#'
#' @examples #
#' # head(eset)[,1:3]
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
  if(dir.exists(dir)) {
    stop(paste0("The directory ", dir, "already exists, do you want to delete the 
         entire directory and re-create it ? ! IF YOU ARE SURE WHAT YOU ARE DOING,
         please use the function mkdir.p\n"))
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
  unlink(dir)
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
#'
#' @return list Large list
#' @export
#' @author Jiang
getGencodeV22.annot <- function() {
  data('gencode.v22.annot')
}