#' @title immuneScore
#' @description Calculate immune score using RNAseq data or mcroarray data
#' 
#' @param exprs matrix, gene expression matrix, NOT log-transformed
#' @param method character, one of "mcp_counter", "quantiseq", "epic", "xcell", "estimate", "timer", "cibersort", "ssGSEA".
#' @param tcga_abbr character, TCGA abbreviation, only used for method = "timer"
#' @param perm number, default 1000(recommend 1000), only used for method = "cibersort"
#' @param QN logical, default F, RNAseq data: recommend QN = F, otherwise T, only used for method = "cibersort"
#' @param absolute absolute logical, if T, return absolute score, only used for method = "cibersort"
#' @param abs_method character, default 'sig.score', can be 'no.sumto1', only used for method = "cibersort"
#' @param gsva_sig_list list, default NULL, only used for method = "ssGSEA"
#'
#' @importFrom immunedeconv deconvolute deconvolute_estimate
#' @importFrom tibble column_to_rownames
#' @importFrom GSVA gsva ssgseaParam
#' @include CIBERSOFT.R
#' 
#' @return a data frame
#' @export
#' 
#' @details
#' The input data is a gene × sample gene expression matrix. \cr
#' In general values should be \cr
#'  TPM-normalized \cr
#'  not log-transformed. \cr
#' For xCell and MCP-counter this is not so important. \cr
#'  xCell works on the ranks of the gene expression only and 
#'  MCP-counter sums up the gene expression values.\cr
#' Rownames are expected to be HGNC gene symbols. \cr
#' Instead of a matrix, immunedeconv also supports ExpressionSets\cr
#' =======================================\cr
#' recommend mcp_counter when comparing between samples
#' EPIC, quanTIseq, CIBERSORT abs : **both**\cr
#' CIBERSORT : **between-cell-type comparisons**\cr
#' MCP-counter, xCell, TIMER, ConsensusTME, ESTIMATE, ABIS: **between-sample comparisons**
#' @examples
#' \dontrun{
#' exprs <- read.table("xx.txt", sep = "\t", row.names = 1, header = T, check.names = F)
#' exprs %>% lzhead
#'  #          TCGA-05-4249-01A TCGA-05-4250-01A TCGA-05-4382-01A
#'  # FGB            84.3223886         17.19138     2495.5995163
#'  # PITX2         127.4491065        125.84088       49.7737300
#'  # SERPINA4        1.2873647       3021.55644        0.0000000
#'  }
immuneScore <- function(exprs, method, tcga_abbr,
                        perm = 1000, QN = F, absolute = T, abs_method='sig.score',
                        gsva_sig_list) {
  if (!is.data.frame(exprs)) {
    stop("expr must be a data.frame")
  }
  if (!is.character(method)) {
    stop("method must be a character")
  }
  deconvolution_methods <- c("mcp_counter", "quantiseq", "epic", "xcell", 
                             "estimate", "timer", "cibersort", "ssGSEA")
  if (!method %in% deconvolution_methods) {
    stop("method must be one of ", paste(deconvolution_methods, collapse = ", "))
  }
  if (method == "mcp_counter") {
    res <- deconvolute(exprs, "mcp_counter")
  } else if (method == "quantiseq") {
    res <- deconvolute(exprs, "quantiseq")
  } else if (method == "epic") {
    res <- deconvolute(exprs, "epic")
  } else if (method == "xcell") {
    res <- deconvolute(exprs, "xcell")
  } else if (method == "estimate") {
    res <- deconvolute_estimate(exprs)
  } else if (method == "timer") {
    res <- deconvolute(exprs, "timer", indications=rep(tcga_abbr, ncol(exprs)))
    res <- res %>% tibble::column_to_rownames("cell_type")
  } else if (method == "cibersort") {
    tmpenv <- new.env()
    data(LM22, package = "LZ", envir = tmpenv)
    LM22 <- tmpenv$LM22
    res <- CIBERSORT(sig_matrix = LM22, mixture_obj = exprs, 
                     perm, QN, absolute, abs_method)
  } else if (method == "ssGSEA") {
    tmpenv <- new.env()
    data(CELL28, package = "LZ", envir = tmpenv)
    CELL28 <- tmpenv$CELL28
    dat <- as.matrix(exprs)
    if(max(na.omit(dat)) > 50) {
      warning("The input data is not log-transformed, the program will log2-transform it.")
      dat <- log2(dat + 1)
    }
    if (missing(gsva_sig_list)) {
      cat("ssGSEA method needs gsva_sig_list, using default CELL28.\n")
      ssgseaPar <- ssgseaParam(dat, CELL28)
      res <- gsva(ssgseaPar)
    } else {
      # 否则请使用内置数据集
      # 翻译：否则请使用内置数据集
      #
      warning("ssGSEA method needs gsva_sig_list, using user-defined.\n")
      warning("Please make sure the gsva_sig_list is appropriate and formatted, otherwise please use the built-in dataset:CELL28\n")
      ssgseaPar <- ssgseaParam(dat, gsva_sig_list)
      res <- gsva(ssgseaPar)
    }
    res <- data.frame(res, check.names = F)
  }
  return(res)
}


#' @title immuneScore_clean
#' @description Clean the immuneScore() result
#' 
#' @param immuneScore_df data frame, the result of immuneScore()
#' @param method character, one of "mcp_counter", "quantiseq", "epic", "xcell", "estimate", "timer", "cibersort", "ssGSEA"
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' res <- immuneScore(res_mcp_counter, method = "mcp_counter")
#' }
immuneScore_clean <- function(immuneScore_df, method) {
  if (method %in% c("mcp_counter", "quantiseq", "epic", "xcell")) {
    df <- immuneScore_df %>% 
      tibble::column_to_rownames(var = names(immuneScore_df)[1]) %>% 
      t() %>% data.frame(check.names = F)
  } else if (method %in% c("estimate", "ssGSEA", "timer")) {
    df <- immuneScore_df %>% t() %>% data.frame(check.names = F)
  } else if (method %in% c("cibersort", "cibersort_abs")) {
    df <- immuneScore_df[, 1:22] %>% data.frame(check.names = F)
  }
  return(df)
}


#' @title immuneScore_onestep
#' @description Calculate immune score using RNAseq data or mcroarray data using multiple methods in one step
#'
#' @param exprs matrix, gene expression matrix, NOT log-transformed
#' @param tcga_abbr character, TCGA abbreviation, only used for method = "timer"
#' @param QN logical, default F, RNAseq data: recommend QN = F, otherwise T, only used for method = "cibersort"
#' @param value logical, default T, if T, assign variables to global environment
#' @param rapid logical, default F, if T, only run mcp_counter, quantiseq, epic, estimate, ssGSEA, xcell, timer\cr
#' otherwise, also run cibersort and cibersort_abs(will take a long time)
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' res_df <- immuneScore_onestep(exprs, tcga_abbr = "brca")
#' }
immuneScore_onestep <- function(exprs, tcga_abbr = "brca", value = T,
                                QN = F, rapid = F) {
  res_suffix <- c("mcp_counter", "quantiseq", "epic", "estimate", "ssGSEA", 
                  "xcell", "timer", "cibersort_abs", "cibersort")
  res_mcp_counter <- immuneScore(exprs = exprs, method = "mcp_counter")
  res_quantiseq <- immuneScore(exprs = exprs, method = "quantiseq")
  res_epic <- immuneScore(exprs = exprs, method = "epic")  #ww
  res_estimate <- immuneScore(exprs = exprs, method = "estimate")
  res_ssGSEA <- immuneScore(exprs = exprs, method = "ssGSEA")
  res_xcell <- immuneScore(exprs = exprs, method = "xcell") #w
  # tcga_abbr为tcga缩写，如brca, coad, luad ...
  res_timer <- immuneScore(exprs = exprs, method = "timer", tcga_abbr = tcga_abbr)
  # 时间长
  if (!rapid) {
    res_cibersort_abs <- immuneScore(exprs = exprs, method = "cibersort", 
                                     perm = 1000, QN = QN, absolute = T, 
                                     abs_method = "sig.score")
    res_cibersort <- immuneScore(exprs = exprs, method = "cibersort", 
                                 perm = 1000, QN = QN, absolute = F, 
                                 abs_method = "sig.score")
    res_suffix <- res_suffix[1:7]
  }
  # assign variables to global environment
  if (value) {
    for (x in res_suffix) {
      assign(paste0("res_", x), get(paste0("res_", x)), envir = .GlobalEnv)
    }
  }
  # integrate
  res_list <- map(res_suffix, ~immuneScore_clean(get(paste0("res_", .x)), .x))
  names(res_list) <- res_suffix
  # 将res.list的每个表格的列名+_他们各自的list的名字
  res_df <- map2_dfc(res_list, names(res_list), ~{
    colnames(.x) <- paste0(.y, "_", colnames(.x))
    return(.x)
  })
  return(res_df)
}