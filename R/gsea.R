#' GSEA plot
#' @description GSEA plot for a result of clusterProfiler::GSEA()
#' 
#' @param GSEA.obj obeject result of  clusterProfiler::GSEA()
#' @param num number the row postion of result of  clusterProfiler::GSEA() 
#'
#' @return plot
#' @export
#' @importFrom enrichplot gseaplot2
#' @importFrom ggsci pal_jco
#'
#' @author Jiang
DEGp_GSEA <- function(GSEA.obj, num=1) {
  p.color <- ifelse(GSEA.obj[,"enrichmentScore"][num] < 0, ggsci::pal_jco()(9)[6],
                    ggsci::pal_jco()(9)[4])
  pic_gsea <- enrichplot::gseaplot2(GSEA.obj, num,
                                    title = GSEA.obj[,"ID"][num],
                                    color= p.color, #线条颜色
                                    base_size = 14, #基础字体的大小
                                    #subplots = 1:2, #展示上2部分
                                    rel_heights = c(1.2, 0.3, 0.5),
                                    pvalue_table = F)
  x.pos <- length(GSEA.obj@geneList) - 2000
  y.pos <- ifelse(GSEA.obj[num,"enrichmentScore"] < 0, -0.1, GSEA.obj[num,"enrichmentScore"]-0.1)
  annotate_p <- paste0("PValue=", round(GSEA.obj[num,"pvalue"], 6), "\nPadjust=", round(GSEA.obj[num,"p.adjust"], 6))
  pic_gsea[[1]] <- pic_gsea[[1]] + annotate("text", x = x.pos, y = y.pos,
                                            label = annotate_p, size = 5)
  return(pic_gsea)
}


#' GSEA plot loop
#' @description Plot all pictures from a result of clusterProfiler::GSEA()
#' 
#' @param gsea.all object result of clusterProfiler::GSEA() using lots pathway gmt files
#' @param result_dir character the outdir name
#' @param xl_filename character the out filename
#'
#' @return NULL
#' @export
#' @importFrom dplyr `%>%`
#' @importFrom dplyr filter
#'
#' @author Jiang
DEGp_GSEA_plotALL <- function(gsea.all, result_dir, xl_filename) {
  result_dir <- dirclean(result_dir)
  mkdir(result_dir)
  for( i in seq_along(gsea.all[,'ID'])) {
    pic_gsea <- DEGp_GSEA(gsea.all, num = i)
    id = gsea.all[,'ID'][i]
    # saveRDS(gsea, file = "result/M1_p0.05fdr0.1log2fc1/gesa.ecm.plot.rds")
    # 使用正则表达式检查字符串是否包含不合法的字符
    if (grepl("[/\\:*?\"|<>]", id)) {
      # 如果包含不合法字符，则进行替换
      st_id <- gsub("[/\\:*?\"|<>]", "_", id)
      pdf(paste0(result_dir, "/gsea_kegg_",st_id, ".pdf"), height = 6, width = 7)
    } else {
      # 如果字符串合法，不进行替换
      pdf(paste0(result_dir, "/gsea_kegg_",id, ".pdf"), height = 6, width = 7)
    }
    print(pic_gsea) # labels = 'A', label_x = 0.05
    dev.off()
    cat(i, id , "\n")
  }
  cat("画图全部完成，正在保存数据中...\n")
  gsea.alldf <- data.frame(gsea.all, check.names = F)
  gsea.alldf.sig <- gsea.alldf %>% filter(., pvalue < 0.05 & p.adjust < 0.1)
  writexl::write_xlsx(list(GSEA.ALL=gsea.alldf, GSEA.sig = gsea.alldf.sig), 
                      path = paste0(result_dir, "/", xl_filename, ".gsea.xlsx"))
  saveRDS(gsea.all, file = paste0(result_dir, "/", xl_filename, ".gsea.rds"))
  cat("保存完成\n")
}


#' Find a pathway gmt
#' @description Find a pathway from the gmt consist of lots of pathway\cr
#' NOTE: this function's result is still a gmt(dataframe)
#' 
#' @param your_pathway character the taget pathway name, support regular expression
#'
#' @return gmt
#' @export
#' @importFrom dplyr `%>%`
#' @importFrom dplyr filter
#'
#' @author Jiang
find_pathway <- function(your_pathway = "^Ferr") {
  if ( sum(grepl(your_pathway, gmt$term)) > 0) {
    pathgmt <- gmt %>% filter(., grepl("^Ferr", term))
    return(pathgmt)
  } else {
    stop("对不起，找不到您指定的通路！")
  }
}


#' GMT Long to Wide data use groub_by
#' @description transformat long GMT dataframe to wide dataframe
#' 
#' @param gmtdf dataframe the long dataframe consist of two coloumns.\cr
#' # term        gene
#' # pathway1   gene1
#' # pahtway1   gene2
#' # ...        ...
#' # pahtwayn   gene1
#' # pathwayx   genex
#'
#' @return gmtdf.w
#' @export
#' @importFrom dplyr `%>%`
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#'
#' @examples
#' # return datafame as this:
#' # term       gene
#' # pathway1   gene1,gene2,...
#' # pathway2   gene1,genex,...
#' @author Jiang
gmt_longTowide2 <- function(gmtdf) {
  gmt.w <- gmtdf %>% 
    dplyr::group_by(term) %>%
    dplyr::summarize(gene = paste(gene, collapse = ", "))
  return(gmt.w)
}


#' GMT Long to Wide data
#' @description trans long GMT dataframe to wide dataframe
#' 
#' @param df dataframe the long dataframe consist of two coloumns.\cr
#' # term        gene
#' # pathway1   gene1
#' # pahtway1   gene2
#' # ...        ...
#' # pahtwayn   gene1
#' # pathwayx   genex
#'
#' @return gmt.wide
#' @importFrom reshape2 dcast
#' @importFrom dplyr `%>%`
#' @importFrom dplyr select
#'
#' @examples
#' # return datafame as this:
#' # term       gene
#' # pathway1   gene1,gene2,...
#' # pathway2   gene1,genex,...
#' @author Jiang
gmt_longTOwide3 <- function(df) {
  dftt <- reshape2::dcast(df, term~gene)
  xn <- names(dftt)[2:ncol(dftt)]
  dftt$GENE <- apply(dftt[, 2:ncol(dftt)], 1, function(x) {
    xn[x>0] %>% sort() %>% paste(., collapse=',')
  })
  df.new <- dftt %>% dplyr::select(term, GENE)
  return(df.new)
}


#' Run GSEA of gmt
#' @description Run GSEA analysis with a gmt file consist of lots of pathway or a single pathway
#' @param genelist list, a name list store the value of all gene after arrage
#' @param gmt_set dataframe or a object, a gmt dataframe filename or gmt object, col1: term, col2: gene
#' #            term  gene
#' # GOBP_PYROPTOSIS  AIM2
#' # GOBP_PYROPTOSIS  APIP
#' # GOBP_PYROPTOSIS CASP1
#' # GOBP_PYROPTOSIS CASP4
#' # ...               ...
#' # KEGG...          ...
#' @param pic.save logical if save the picture, default is F, it's not reccommed set as T if gmt is lots of pathway
#' @param outdir character the picture outdir
#' @param filename character the picture filename
#' @param pic.w number picture width
#' @param pic.h number picture height
#' @param num number the row postion of result of  clusterProfiler::GSEA() 
#' @param pvalue number the p value of gsea analysis
#'
#' @return gsea
#' @export
#' @importFrom clusterProfiler GSEA
#'
#' @author Jiang
DEG_runGSEA <- function(genelist, gmt_set, pic.save, outdir, filename,
                    pic.w=7, pic.h=6, num = 1, pvalue = 1) {
  stopifnot(is.character(gmt_set) | is.data.frame(gmt_set))
  if (is.character(gmt_set)) {
    gmt <- read.gmt(gmt_set) 
  } else {
    gmt <- gmt_set
  }
  #
  gsea <- clusterProfiler::GSEA(genelist, TERM2GENE = gmt, pvalueCutoff = pvalue) #GSEA分析
  if (pic.save) {
    stopifnot(!missing(outdir) | !missing(filename))
    pic_gsea <- DEGp_GSEA(gsea.all, num = num)
    outdir <- dirclean(outdir)
    mkdir(outdir)
    pdf(paste0(outdir, "/", filename, ".pdf"), width = pic.w, height = pic.h)
    print(pic_gsea)
    dev.off()
  }
  return(gsea)
}