#' GO analysis
#' @description GO analysis by ENTREZID, plot go tree picture and write to disk, retern the result
#'
#' @param genelist dataframe, the dataframe must contain "ENTREZID" coloum
#' @param orgdb character, org type
#' @param sigNodes number, the tree plot nodes
#' @param resultdir character, the output dir
#' @param filemark character, special part of output filename
#'
#' @return # list, GO analsis object go_list
#' @export
#' @importFrom clusterProfiler plotGOgraph enrichGO
#'
#' @examples #
DEG_GO <- function(genelist, orgdb="org.Hs.eg.db", sigNodes=20,
                   resultdir, filemark) {
  ## go分析，并画go关系图，保存golist对象
  #genelist为数据框，需有一列为entrzid
  go_list <- list()
  for (ont in c('ALL', 'BP', 'CC', 'MF')) {
    GO <- enrichGO(genelist[, "ENTREZID"],#GO富集分析BP模块
                   OrgDb = orgdb,
                   keyType = "ENTREZID",
                   ont = ont,
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 1,
                   readable = T)
    if (ont=='ALL') {
      go_list[[ont]] <- GO
      go_all_df <- data.frame(GO)
      #write.table(go_all_df, sep = '\t', quote = F,
      #            file = paste0(resultdir, '/GOall', '_',filemark, '.txt'))
    } else {
      go_list[[ont]] <- GO
      # plot关系图
      fname <- paste0(resultdir, '/', ont, '_',filemark, '.pdf')
      cat(fname)
      pdf(fname, height = 8, width = 12)
      plotGOgraph(GO, firstSigNodes = sigNodes)
      dev.off()
    }
  }
  # save(go_list, file = paste0(resultdir, '/GOall', '_',filemark, '.Rdata'), compress = T)
  return(go_list)
}


#' KEGG analysis
#' @description KEGG analysis by UNIPROT, retern the result
#'
#' @param genelist dataframe, the dataframe must contain "UNIPROT" coloum
#' @param orgdb character, org type
#' @param org_kegg character, kegg org type
#'
#' @return # list, KEGG analsis object & its df, list(KEGG=KEGG, pSigDF=keggdf)
#' @export
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom DOSE setReadable
#' @importFrom dplyr filter
#'
#' @examples #
DEG_KEGG <- function(genelist, orgdb="org.Hs.eg.db", org_kegg='hsa') {
  # resultdir, filemark) {
  ## 自定义函数，kegg分析，保存kegg对象（原y叔函数实为p.adjust）,人为取出p<0.05的数据
  #genelist为数据框，需有一列为entrzid
  KEGG <- enrichKEGG(genelist[, "UNIPROT"], keyType = 'uniprot',
                     organism = org_kegg, minGSSize = 1,
                     pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG <- setReadable(KEGG, OrgDb = orgdb, keyType="UNIPROT")
  keggdf <- filter(KEGG@'result', pvalue < 1)
  #save(KEGG, file = paste0(resultdir, '/KEGGpathway', '_', filemark, '.Rdata'), compress = T)
  #write.table(keggdf, sep = '\t', quote = F,
  #           file = paste0(resultdir, '/KEGG', '_', filemark, '.txt'))
  return(list(KEGG=KEGG, pSigDF=keggdf))
}


#' GO analysis data prepare
#' @description Transform the DEG anlysis object(all_father$DIFF.ALL') to GO analsis
#' Note: The go enrich result of the length of different deg genes is different!
#'
#' @param resdf dataframe, the DEG anlysis object(all_father$DIFF.ALL')
#' @param logfc number, logfc
#' @param p number, pvale
#' @param fdr number, fdr
#'
#' @return # list(all=all, up=up, down=down)
#' @export
#'
#' @examples #
DEG_prepareGOglist <- function(resdf, logfc, p=0.05, fdr=0.1) {
  up <- resdf %>% filter(log2FC > logfc & PValue < p & FDR < fdr) %>% .[,1]
  down <- resdf %>% filter(log2FC < -logfc & PValue < p & FDR < fdr) %>% .[,1]
  all <- c(up, down)
  return(list(all=all, up=up, down=down))
}


#' run GO analysis in batches(upgene, downgene, allgene)
#' @description run GO analysis in batches(upgene, downgene, allgene) and write result to xlsx
#'
#' @param outdir character, the output dir
#' @param genelist list, go genelist, result of DEG_prepareGOglist
#'
#' @return # GO object & write result to xlsx
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom xlsx write.xlsx
#'
#' @examples #
DEG_runGO <- function(outdir, genelist) {
  # 对logFC迭代，每次FC新建一个目录，下存upgene, downgene, allgene的GO结果
  mapply(function(x, y) {
    outd <- paste0(outdir, "/FC_", y)
    mkdir(outd)

    mapply(function(xx, yy) {
      suppressMessages({ suppressWarnings({
          gene_df <- bitr(xx, fromType = 'SYMBOL', toType = 'ENTREZID',
                          OrgDb = GO_database)
        }) })

      f_mark = paste0(y, "_", yy)
      cat(f_mark, "\n", "......", "\n")
      go.an <- DEG_GO(genelist = gene_df, resultdir = outd, filemark = yy)
      write.xlsx(go.an$ALL, file = paste0(outd, '/go.', y,'.xlsx'), col.names = T,
                 row.names = F, sheetName=yy, append = T)
    }, x, names(x))

    cat("  Done!\n")
  }, genelist, names(genelist), SIMPLIFY = F)
}

#' run KEGG analysis in batches(upgene, downgene, allgene)
#' @description run KEGG analysis in batches(upgene, downgene, allgene) and write result to xlsx
#'
#' @param outdir character, the output dir
#' @param genelist list, go genelist, result of DEG_prepareGOglist
#'
#' @return # KEGG object & write result to xlsx
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom xlsx write.xlsx
#'
#' @examples #
DEG_runKEGG <- function(outdir, genelist) {
  # 对logFC迭代，每次FC新建一个目录，下存upgene, downgene, allgene的KEGG结果
  mapply(function(x, y) {
    outd <- paste0(outdir, "/FC_", y)
    mkdir(outd)

    mapply(function(xx, yy) {
      suppressMessages({ suppressWarnings({
          gene_df <- bitr(xx, fromType = 'SYMBOL', toType = 'UNIPROT',
                          OrgDb = GO_database)
        }) })

      f_mark = paste0(y, "_", yy)
      cat(f_mark, "\n", "......", "\n")
      write.xlsx(gene_df, file = paste0(outd, '/genelist.', y,'.xlsx'),
                 col.names = T, row.names = F, sheetName=yy, append = T)
      kegg.an <- DEG_KEGG(genelist = gene_df)
      write.xlsx(kegg.an$pSigDF, file = paste0(outd, '/kegg.', y,'.xlsx'),
                 col.names = T, row.names = F, sheetName=yy, append = T)
    }, x, names(x))

    cat("  Done!\n")
  }, genelist, names(genelist), SIMPLIFY = F)
}


