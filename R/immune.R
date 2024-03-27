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