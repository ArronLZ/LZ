#' GO analysis
#' @description GO analysis by ENTREZID, plot go tree picture and write to disk, retern the result
#'
#' @param genelist dataframe, the dataframe must contain "ENTREZID" coloum
#' @param orgdb character, org type
#' @param sigNodes number, the tree plot nodes
#' @param resultdir character, the output dir
#' @param filemark character, special part of output filename
#' @param rapid logical, default is T, meaning not plot tree plot
#'
#' @return # list, GO analsis object go_list
#' @export
#' @importFrom clusterProfiler plotGOgraph enrichGO
#'
#' @examples #
#' @author Jiang
DEG_GO <- function(genelist, orgdb="org.Hs.eg.db", sigNodes=20,
                   resultdir, filemark, rapid = T) {
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
      #go_all_df <- data.frame(GO)
    } else {
      go_list[[ont]] <- GO
      # plot关系图
      if (!rapid) {
        resultdir <- dirclean(resultdir)
        mkdir.p(resultdir)
        fname <- paste0(resultdir, '/', ont, '_', filemark, '.pdf')
        cat(fname)
        pdf(fname, height = 8, width = 12)
        plotGOgraph(GO, firstSigNodes = sigNodes)
        dev.off()
      }
    }
  }
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
#' @author Jiang
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
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' 
#' @author Jiang
DEG_prepareGOglist <- function(resdf, logfc, p=0.05, fdr=0.1) {
  up <- resdf %>% filter(log2FC > logfc & PValue < p & FDR < fdr) %$% Gene
  down <- resdf %>% filter(log2FC < -logfc & PValue < p & FDR < fdr) %$% Gene
  all <- c(up, down)
  return(list(all=all, up=up, down=down))
}


#' Enrich analysis data prepare
#' @description Transform the DEG anlysis object(all_father$DIFF.ALL') to GO analsis
#' 
#' @param resdf dataframe, the DEG anlysis object(all_father$DIFF.ALL')
#' @param fc.list list, name list such as list('1.5' = log2(1.5), '2' = log2(2))
#' @param p number, pvale
#' @param fdr number, fdr
#'
#' @return list
#' @export
#'
#' @author Jiang
DEG_prepareENRICH <- function(resdf, fc.list, p = 0.05, fdr = 0.1) {
  gogenelist <- lapply(fc.list, function(x) {
    DEG_prepareGOglist(resdf, logfc = x, p = p, fdr = fdr) })
  names(gogenelist) <- paste0("p", p, ";q", fdr, ";fc", names(gogenelist))
  return(gogenelist)
}


#' run KEGG analysis inFCs(2, 1.5) in batches(upgene, downgene, allgene)
#' @description run GO analysis in batches(upgene, downgene, allgene) and write result to xlsx
#'
#' @param genelist list, must be two-dimensional name list, result of sapply DEG_prepareGOglist
#' @param outdir character, the output dir
#' @param glist.save logical default is TRUE, wether save genelist
#' @param rungo logical default is TRUE, wether run go analysis
#' @param runkegg logical default is TRUE, wether run kegg analysis
#' @param rapid loggicla defaut is F, if set F, the function will plot go tree picture, which make cost much time.
#' @param orgdb character default is "org.Hs.eg.db"
#' @param org_kegg character default is 'hsa'
#' @param sigNodes number default is 20
#'
#' @return # GO object & write result to xlsx
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom writexl write_xlsx
#'
#' @examples #
#' @author Jiang
DEG_runENRICH <- function(genelist, outdir, glist.save = T, 
                          rungo=T, runkegg=T, rapid = T,
                          orgdb = "org.Hs.eg.db", org_kegg='hsa',
                          sigNodes=20
                          ) {
  # 对logFC迭代，每次FC新建一个目录，下存upgene, downgene, allgene的KEGG结果
  kegg.list <- list()
  kegg.df.list <- list()
  #
  go.list <- list()
  go.df.list <- list()
  
  name.list <- names(genelist)
  progressbar <- (1:length(name.list))/length(name.list)
  progressbar <- paste0(round(progressbar*100) - 0.01, "%")
  i = 1
  cat("========", "START ENRICH...", "\n")
  for (name_x in name.list) {
    cat("  ", name_x, "...\n")
    cat("========\n")
    genedf.list <- sapply(genelist[[name_x]], function(x) {
      clusterProfiler::bitr(x, fromType = 'SYMBOL', 
                            toType = c('UNIPROT', 'ENTREZID', 'ENSEMBL'),
                            OrgDb = orgdb) 
    }, simplify = F)
    #
    outd <- paste0(dirclean(outdir), "/", name_x)
    mkdir.p(outd)
    #
    if (glist.save) {
      genedf.list.symbol <- sapply(genedf.list, function(x) {
        x = x[!duplicated(x$SYMBOL),]
        x[, 1, drop = F]
      }, simplify = F)
      writexl::write_xlsx(genedf.list, 
                          path = paste0(outd, "/genelist.multiID.", name_x, ".xlsx"))
      writexl::write_xlsx(genedf.list.symbol, 
                          path = paste0(outd, "/genelist.symbol.", name_x, ".xlsx"))
    }
    if (runkegg) {
      # part KEGG
      cat("  ... KEGG START...\n")
      # genelist, orgdb="org.Hs.eg.db", org_kegg='hsa'
      kegg.an <- sapply(genedf.list, function(x) 
        DEG_KEGG(genelist = x, orgdb=orgdb, org_kegg=org_kegg), simplify = F)
      # 向kegg.list保存kegg分析结果
      kegg.list[[name_x]] <- kegg.an
      # 获取kegg分析结果中的dataframe
      kegg.an.df <- sapply(kegg.an, function(x) x$pSigDF, simplify = F)
      # 向kegg.df.list保存kegg dataframe结果
      kegg.df.list[[name_x]] <- kegg.an.df
      writexl::write_xlsx(kegg.an.df, 
                          path = paste0(outd, "/kegg.", name_x, ".xlsx"))
      cat("  ... KEGG END!\n")
    }
    # part GO 
    if (rungo) {
      cat("  ... GO START...\n")
      go.an <- mapply(function(x, xn) {
        cat("  ... ", xn, "...\n")
        DEG_GO(genelist = x, orgdb=orgdb, sigNodes=sigNodes, 
               resultdir = outd, filemark = paste0(name_x, "_", xn), 
               rapid = rapid)
      }, genedf.list, names(genedf.list), SIMPLIFY = F)
      go.list[[name_x]] <- go.an
      go.an.df <- sapply(go.an, function(x) data.frame(x$ALL), simplify = F)
      go.df.list[[name_x]] <- go.an.df
      writexl::write_xlsx(go.an.df, 
                          path = paste0(outd, "/go.", name_x, ".xlsx"))
      cat("  ... GO END!\n")
    }
    cat("==============", progressbar[i], "\n")
    i = i + 1
  }
  #
  return.list <- list('KEGG'= kegg.list, 'KEGGDF'=kegg.df.list,
                      'GO'=go.list, 'GODF'=go.df.list)
  cat(" ALL Done! 100% \n")
  return(return.list)
}
