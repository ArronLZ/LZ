#' Gene co-expression analysis
#' @description This function is used to calculate the co-expression of the target gene \cr
#' with all other genes in the expression matrix.
#'
#' @param expr data.frame, the expression matrix dataframe
#' @param taget_gene character, the gene name to be used as the target gene
#'
#' @return data.frame, the co-expression result
#' @export
#' @import dplyr
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \dontrun{
#' expr %>% lzhead()
#' === shape ===
#'   19057 rows 113 columns
#' === head 1 : 4 ===
#'          TCGA-A1-A0SK TCGA-A1-A0SP TCGA-A2-A04U TCGA-A2-A0CM
#' TSPAN6      12.647009    12.064743    12.216746    10.954196
#' TNMD         0.000000     2.807355     4.392317     1.584963
#' FGR          6.882643     9.856426     7.169925     9.668885
#' res_cor <- gene_co_expression(expr, "PHLDB2")
#' }
gene_co_expression <- function(expr, taget_gene = "PHLDB2") {
  res <- apply(expr, 1, function(x) {
    cor.res <- cor.test(x, as.numeric(expr[taget_gene,]))
    c(cor.res$p.value, cor.res$estimate)
  }) %>% t() %>% data.frame() %>% 
    setNames(c("PValue", "Cor")) %>% 
    mutate(FDR = p.adjust(PValue, method = "fdr")) %>% 
    tibble::rownames_to_column("Gene")
  return(res)
}
