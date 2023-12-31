#' process xena data(annot & remove duplicate &or fpkm to tpm)
#' @description process xena data: 1. annot symbol 2. removes duplicate row(keep the max rowmean symbol) 3. fpkm to tpm(optional)
#' @param eset dataframe refer to examples
#' @param annot dataframe refer to examples, must contain the coloumn "id" and "gene"
#' @param log logical default is TRUE, if the eset have been log transformat
#' @param tpm logical when set to TURE eset must be FPKM data
#' @param oid character optional, the xena data first column name usually be Ensembl_ID, if the first column of your data is not Ensemble_ID please modify the param.
#' @param id.simplify logical, default is FALSE, if the id is full ensemble id that have the . , the param show set be TRUE
#' 
#' @return dataframe(1st col is symbol, 2nd col is the ensemble_id), refer to examples
#' @export
#' @import dplyr
#' @importFrom stringr str_split
#'
#' @examples
#' # eset
#' #          Ensembl_ID  TCGA-97-7938-01A  TCGA-55-7574-01A
#' #  ENSG00000000003.13          10.98939          9.967226
#' #   ENSG00000000005.5           4.00000          0.000000
#' #  ENSG00000000419.11          10.25385          9.541097
#' # annot
#' #                 id       gene  chrom  chromStart
#' #  ENSG00000223972.5    DDX11L1   chr1       11869
#' #  ENSG00000227232.5     WASH7P   chr1       14404
#' #  ENSG00000278267.1  MIR6859-3   chr1       17369
#' # return
#' #     SYMBOL        Ensembl_ID  TCGA-MP-A4T4-01A
#' #    5S_rRNA   ENSG00000271924           1.89573
#' #  5_8S_rRNA   ENSG00000275877           0.00000
#' #        7SK   ENSG00000271394           0.00000
XENApre_annot <- function(eset, annot, log = TRUE, tpm=FALSE, 
                       oid="Ensembl_ID", id.simplify = FALSE) {
  if (log == TRUE) {
    eset[,2:ncol(eset)] <- round(2^eset[,2:ncol(eset)] - 1)
  }
  if (tpm == TRUE) {
    cat("提示: 当设置tpm = TURE时必须保证eset数据为FPKM数据！\n")
    # eset[,2:ncol(eset)] <- apply(eset[,2:ncol(eset)], 2, function(x) round(2^x - 1, 0))
    eset[,2:ncol(eset)] <- apply(eset[,2:ncol(eset)], 2, fpkmToTpm)
    #eset[,2:ncol(eset)] <- apply(eset[,2:ncol(eset)], 2, function(x) log2(x + 1))
  }
  if (id.simplify == TRUE) {
    # 简化ensambleID
    eset[, oid] <- str_split(eset[, oid], "\\.", simplify = T)[,1]
  }
  # 注释
  eset <- cbind.data.frame(
    SYMBOL = annot[match(eset[, oid], annot$id), "gene"],
    eset)
  # 去重
  eset$MEAN <- abs(rowMeans(eset[,3:ncol(eset)]))
  eset <- eset %>% arrange(SYMBOL, desc(MEAN))
  eset <- eset[!duplicated(eset$SYMBOL), ]
  rownames(eset) <- NULL
  if (!is.null(eset$MEAN)) { eset$MEAN <- NULL }
  return(eset)
}





