#' Combine method
#' @description combine coloumn from dataframe or colnames
#' 
#' @param col.1 character(n) or numerber(n)
#' @param col.2 character(n) or numerber(n)
#' @param sample1 number the sample number of group1
#' @param sample2 number the sample number of group1
#' @param eset data.frame the expr dataframe
#'
#' @return dataframe or matrix
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom utils combn
#' @importFrom stringr str_split
#'
#' @author Jiang
combn_sample <- function(col.1, col.2, sample1 = 2, sample2 =2, eset) {
  stopifnot(is.numeric(col.1) | is.character(col.1), 
            is.numeric(col.2) | is.character(col.2))
  if (!missing(eset)) {stopifnot(is.data.frame(eset))}
  if (!missing(eset) & is.character(col.1) & is.character(col.2)) {
    stop("若设置了eset参数，则col.1, col.2参数只能为数值型，请改为对应的列位置数值")
  }
  if (missing(eset)) {
    a = combn(col.1, sample1)
    b = combn(col.2, sample2)
  } else {
    a = combn(colnames(eset)[col.1], sample1)
    b = combn(colnames(eset)[col.2], sample2)
  }
  
  aa <- apply(a, 2, function(x) {
    paste(x, collapse = ",")
  })
  bb <- apply(b, 2, function(x) {
    paste(x, collapse = ",")
  })
  com_df <- expand.grid(aa, bb)
  com_df <- cbind(stringr::str_split(com_df$Var1, ",", simplify = T),
                  stringr::str_split(com_df$Var2, ",", simplify = T))
  if (!missing(eset)) {
    com_df <- apply(com_df, 1, function(x) {
      eset %>% dplyr::select(any_of(x))
    })
  }
  return(com_df)
}