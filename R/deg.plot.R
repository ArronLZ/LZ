#' Valcano with mask gene
#' @description plot valcano with mask gene
#'
#' @param result dataframe, all info of DEG analysis
#' @param logFC logfc, default 1
#' @param adj_P padj, default 0.2
#' @param label_geneset character, the mask gene list
#'
#' @return #
#' @export
#' @import ggrepel
#' @import ggsci
#' @import dplyr
#'
#' @examples #
DEGp_Volcano <- function(result, logFC = 1, adj_P = 0.2, label_geneset = NULL) {
  jco <- ggsci::pal_npg()(9)
  result$Type = "NONE"
  result$Type[which(result$log2FoldChange > logFC & result$padj < adj_P)] = "UP"
  result$Type[which(result$log2FoldChange < (-logFC) & result$padj < adj_P)] = "DOWN"
  xlim = max(abs(result$log2FoldChange))
  if(is.null(label_geneset)) {
    p = ggplot(result, aes(x = log2FoldChange, y = -log10(padj)))+
      geom_point(data = result, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
      theme_bw() +
      geom_vline(xintercept = c(-logFC,logFC), lty = 2)+
      geom_hline(yintercept = c(-log10(adj_P)), lty = 2)+
      scale_x_continuous(limits = c(-xlim, xlim))+
      coord_fixed(ratio = ( 2*xlim )/(max(-log10(result$padj[result$padj != 0]), na.rm = T)) )+
      theme(panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"))+
      xlab("log2(Fold Change)")+ylab("-log10(Pvalue)")+
      scale_color_manual(values = c("NONE" = "grey","UP" = "red","DOWN" = "blue"))
  } else {
    p = ggplot(result,aes(x = log2FoldChange, y = -log10(padj)))+
      geom_point(data = result,
                 aes(x = log2FoldChange, y = -log10(padj), color = Type), alpha=0.9)+
      geom_point(data = result[which(result$Row.names %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj)),color = "black",size = 4)+
      geom_point(data = result[which(result$Row.names %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj)),color = "white",size = 2.5)+
      geom_point(data = result[which(result$Row.names %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj),color = Type),size = 1.5)+
      geom_text_repel(data = result[which(result$Row.names %in% label_geneset),],
                      aes(x = log2FoldChange, y = -log10(padj), label = Row.names),
                      fontface = "bold", max.overlaps = 30) + #fontface = "italic"
      theme_bw() +
      geom_vline(xintercept = c(-logFC,logFC),lty = 2)+
      geom_hline(yintercept = c(-log10(adj_P)),lty = 2)+
      scale_x_continuous(limits = c(-xlim,xlim))+
      coord_fixed(ratio = ( 2*xlim )/(max(-log10(result$padj[result$padj != 0]), na.rm = T))  )+
      theme(panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"))+
      xlab("log2(Fold Change)")+ylab("-log10(Pvalue)")+
      scale_color_manual(values = c("UP" = jco[8],"DOWN" = jco[4],"NONE" = "grey"))
  }
  return(p)
}


#' Heatmap
#' @description plot heatmap with selected gene
#'
#' @param df dataframe or matrix, the exprs data(rowname is gene name)
#' @param gene character, the gene name that heatmap show
#' @param pheno dataframe, the group info
#' @param outdir character, dir name
#' @param f_mark character, special part of filename
#' @param nlevel number, the group's factor level
#' @param pic_w number, picture width
#' @param pic_h number,  picture height
#' @param c_w number,  ceil width
#' @param c_h number, ceil height
#' @param angle character, col angle
#' @param cluster_rows logical, default is T
#' @param c_colors character(n), default isc ("blue", "#FFFFCC", "red")
#' @param color.num number, default is 100
#' @param fontsize number, row(sample) annot font size
#' @param fontsize.col number, col(gene) annot font size
#' @param clustering_method character, one of c("ward.D", "ward.D2", "single", "complete", 
#' "average", "mcquitty", "median", "centroid")
#' @param clustering_distance_rows character, one of c("correlation", "euclidean", "maximum", 
#' "manhattan", "canberra", "binary", "minkowski")
#' @param clustering_distance_cols character, as the same to `clustering_distance_rows`
#' @param cluster_cols logical, default is F
#'
#' @return #
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom ggsci pal_d3
#' @import RColorBrewer
#'
#' @examples
#' \dontrun{
#'   groupdf <- data.frame(row.names = colnames(sig.expr),
#'                         Type = c(rep("11A", 3), rep("01A", 3)),
#'                         check.names = F)
#'   groupdf[,1] <- as.factor(groupdf[,1])
#'   m <- c("ward.D", "ward.D2", "single", "complete", 
#'           "average", "mcquitty", "median", "centroid")
#'   d <- c("correlation", "euclidean", "maximum", "manhattan",
#'          "canberra", "binary", "minkowski")
#'   for (i in m) {
#'     for (i.d in d) {
#'       cat(i, i.d, "\n")
#'       DEGp_Pheat(
#'         sig.expr,
#'         gene = show.gene, 
#'         pheno = groupdf,
#'         outdir = "D:/xxx/yyy",
#'         f_mark = paste0("heatmap.", i, "_", i.d),
#'         nlevel = 2,
#'         pic_w = 5,
#'         pic_h = 8,
#'         c_w = 18,
#'         c_h = 10,
#'         fontsize = 10,
#'         fontsize.col = 14,
#'         angle = "90",
#'         cluster_rows = T,
#'         cluster_cols = T,
#'         clustering_method = i,
#'         clustering_distance_rows=i.d
#'       )  
#'     }
#'   }
#' }
DEGp_Pheat <- function(df, gene, pheno, outdir, f_mark = "", nlevel = 8, pic_w = 6, 
                       pic_h = 8, c_w = 6, c_h = 10, c_colors = c("blue", "#FFFFCC", "red"),
                       color.num = 100, fontsize = 12, fontsize.col = 8, angle = "90", 
                       cluster_rows = T, cluster_cols = F, clustering_method = "complete",
                       clustering_distance_rows="euclidean", clustering_distance_cols="euclidean") {
  heat_matrix <- df[rownames(df) %in% gene, ,]
  matdf <- pheno
  if (nlevel != 2) {
    mat_colors <- list((ggsci::pal_d3())(nlevel))
  }
  else {
    mat_colors <- list(c("dimgrey", "darkorange"))
  }
  names(mat_colors)[1] <- names(matdf)[1]
  names(mat_colors[[1]]) <- levels(matdf[, 1])
  col <- colorRampPalette(brewer.pal(12, "Paired"))(24)
  mkdir(dirclean(outdir))
  # dev.off()
  fname <- paste0(dirclean(outdir), "/heatmap.", f_mark, ".pdf")
  pdf(fname, width = pic_w, height = pic_h)
  p <- pheatmap(heat_matrix, scale = "row", 
                clustering_method = clustering_method,
                clustering_distance_rows = clustering_distance_rows,
                clustering_distance_cols = clustering_distance_cols,
                cluster_rows = cluster_rows, 
                cluster_cols = cluster_cols, 
                color = colorRampPalette(colors = c_colors)(color.num), 
                annotation_col = matdf, annotation_colors = mat_colors, 
                show_colnames = T, show_rownames = T, border_color = NA, 
                cellwidth = c_w, cellheight = c_h, fontsize = fontsize, 
                fontsize_col = fontsize.col, 
                angle_col = angle, drop_levels = T)
  print(p)
  dev.off()
  #png(filename = paste0(outdir, "/heatmap.", f_mark, ".png"), 
  #    height = pic_h, width = pic_w, units = "in", res = 1000)
  #print(p)
  #dev.off()
}


#' Venn Plot
#' @description Plot Venn 2-5 group data
#'
#' @param datalist list, such as list(diffan$deg$symbol, diffagek, ...), no more than 5-dimensional-list
#' @param transparency number, transparency, 0-1, also can be NULL, transparency = NULL means no fill color
#' @importFrom ggsci pal_jco
#' @importFrom scales alpha
#' @importFrom VennDiagram venn.diagram
#' 
#' @return plot
#' @export
#'
#' @examples
#' \dontrun{
#' lzpv <- list(a=1:10, b=5:16, c=15:23, d=c(1,2,4,8,8,16:19,30:34))
#' p <- lzplot_veen(datalist = lzpv, transparency = 0.6)
#' }
#' @author Jiang
DEGp_venn <- function(datalist, transparency = 0.6) {
  # datalist为list,见例子 。 注意本函数返回结果必须使用grid.draw()函数打印
  # transparency 为 0-1，表示透明度，或设置为NULL，表示不绘制区域颜色
  # 例子：
  # lzpv <- list(a=1:10, b=5:16, c=15:23, d=c(1,2,4,8,8,16:19,30:34))
  # p <- lzplot_veen(datalist = lzpv)
  # # # 注意要保存为pdf必须使用 grid.draw()函数打印，长宽比注意控制
  # pdf(file = "./02.diff.geneset/r2.venn_DEGAGES3.pdf", width = 7, height = 6)
  # grid.draw(p)
  # dev.off()
  jco <- ggsci::pal_jco()(9) #jco %>% scales::show_col()
  t <- transparency
  if (length(datalist) == 2) {
    myCol <- c(jco[1], jco[4])
    myCols <- c(alpha(jco[1], t), alpha(jco[4], t))
  } else if (length(datalist) == 3) {
    myCol <- c(jco[1], jco[4], jco[9])
    myCols <- c(alpha(jco[1], t), alpha(jco[4], t), alpha(jco[9], t))
  } else if (length(datalist) == 4) {
    myCol <- c(jco[9], jco[1], "orange", jco[6])
    myCols <- c(alpha(jco[9], t), alpha(jco[1], t),
                alpha("orange", t), alpha(jco[6], t))
  } else if (length(datalist) == 5) {
    myCol <- c(jco[1], jco[2], jco[3], jco[4], jco[9])
    myCols <- c(alpha(jco[1], t), alpha(jco[2], t), alpha(jco[3], t),
                alpha(jco[4], t), alpha(jco[9], t))
  } else {
    stop("the datalist length must be in the range 2 to 5")
  }
  #
  if (is.null(transparency)) { myCols = NULL }
  #
  if ( is.null(names(datalist)) ) {
    category.n = rep('', length(datalist))
  } else {
    category.n = names(datalist)
  }

  grid.newpage()
  p <- venn.diagram(
    x = datalist,
    category.names = category.n, #c("DEGs", "AGE GENE"),
    filename = NULL, #'./02.diff.geneset/r2.venn_DEGAGES.tiff',
    output=TRUE, scale =F, resolution = 2500,
    height = 6000, width = 6000,
    disable.logging = T,
    lwd = 1,
    col = "transparent",
    fill = myCols,
    cex = 1,
    fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(12, -10),
    cat.col = myCol
  )
  grid.draw(p)
  return(p)
}
