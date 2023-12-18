#' prepare the data that valcano need
#' @description prepare the data that valcano need
#'
#' @param df_valcano dataframe, the all data of DEG analysis
#' @param filterc character, the way of show picture
#'
#' @return #
#' @export
#'
#' @examples #
DEGp_prepareVolcano <- function(df_valcano, filterc = 'fppadj') {
  df_valcano <- df_valcano %>%
    dplyr::select(Gene, log2FC, PValue, FDR) %>%
    dplyr::rename(Row.names=Gene,
                  log2FoldChange=log2FC,
                  pvalue = PValue,
                  padj=FDR) %>%
    na.omit()
  if (filterc == 'fppadj') {
    df_valcano$sign <- ifelse(df_valcano$padj < ffdr & df_valcano$pvalue < fpval,
                              ifelse(abs(df_valcano$log2FoldChange) > flogfc,
                                     ifelse(df_valcano$log2FoldChange > flogfc, 'Up', 'Down'),
                                     'Not'),
                              'Not')
  } else if (filterc == 'fpadj') {
    df_valcano$sign <- ifelse(df_valcano$padj < ffdr,
                              ifelse(abs(df_valcano$log2FoldChange) > flogfc,
                                     ifelse(df_valcano$log2FoldChange > flogfc, 'Up', 'Down'),
                                     'Not'),
                              'Not')
  } else {
    df_valcano$sign <- ifelse(df_valcano$pvalue < fpval,
                              ifelse(abs(df_valcano$log2FoldChange) > flogfc,
                                     ifelse(df_valcano$log2FoldChange > flogfc, 'Up', 'Down'),
                                     'Not'),
                              'Not')
  }
  return(df_valcano)
}


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


#' Valcano2 with mask gene
#' @description plot valcano with mask gene
#' @param resdf dataframe the all data of DEG analysis
#' @param filterc character the way of show picture
#' @param pvalue number pvalue, default 0.05
#' @param fdr number fdr, default 0.1
#' @param logfc.p number log2fc positive, default 1
#' @param logfc.n number log2fc negative, default -1
#' @param label_geneset character(n) the genelist need be shown in the valcano, default NULL
#' @param outdir character the picture outdir, default "result"
#' @param filename.base character the picture filename base, default "DEG_xx"
#' @param pic.w number the picture width, default 7
#' @param pic.h number the picture height, default 7
#' @param color.not character the picture color to not sig gene, default "#BEBEBE"
#' @param color.up character the picture color to up sig gene, default "#AC001D"
#' @param color.down character the picture color to down sig gene, default "#3D71AD"
#' @param ymax number the max ylim, default NULL
#' @param size number the picture point size, default 1.5 ...
#' @param shape number the picture point shape, default 16, or 1 ...
#' @param alpha number the transpate of picture, default 0.9
#'
#' @return list
#' @export
#' @import ggrepel
#' @import ggsci
#' @import dplyr
#' @include tool.R
#' @author Jiang
DEGp_Volcano2 <- function(resdf, filterc = "both", 
                          pvalue=0.05, fdr=0.1, logfc.p=1, logfc.n=-1, 
                          label_geneset = NULL, outdir="result", 
                          filename.base = "DEG_xx", pic.w = 7, pic.h =7,
                          color.not = "#BEBEBE", color.up = "#AC001D", 
                          color.down = "#3D71AD", ymax = NULL, size=1.5, 
                          shape=16, alpha = 0.9) {
  mkdir(dirclean(outdir))
  if (!all(c("Gene", "log2FC", "PValue", "FDR") %in% names(resdf) )) {
    stop("数据表格必须包含Gene, log2FC, PValue, FDR这四列")
  }
  if (max(resdf$log2FC) > 30) {
    cat("警告，严重警告：基因的变化倍数列数据最大值大于30,意味着该数据可能尚未log2化\n")
  }
  df <- resdf %>% dplyr::select(Gene, log2FC, PValue, FDR)
  
  if (filterc == "both") {
    df$Category <- ifelse(df$FDR < fdr & df$PValue < pvalue & (
      df$log2FC > logfc.p | df$log2FC < logfc.n), 
      ifelse(df$log2FC > logfc.p, "UP", "DOWN"), "NOT")
  } else if (filterc == "fdr") {
    df$Category <- ifelse(df$FDR < fdr & (df$log2FC > logfc.p | df$log2FC < logfc.n), 
                          ifelse(df$log2FC > logfc.p, "UP", "DOWN"), "NOT")
  } else if (filterc == "pvalue") {
    df$Category <- ifelse(df$PValue < pvalue & (df$log2FC > logfc.p | df$log2FC < logfc.n), 
                          ifelse(df$log2FC > logfc.p, "UP", "DOWN"), "NOT")
  } else { # 为了兼容老版本而设置的，即filterc设置错误的话会按照both的方法进行
    df$Category <- ifelse(df$FDR < fdr & df$PValue < pvalue & (
      df$log2FC > logfc.p | df$log2FC < logfc.n), 
      ifelse(df$log2FC > logfc.p, "UP", "DOWN"), "NOT")
  }
  
  # return(df_valcano)
  xlim = max(abs(df$log2FC))
  #minp <- -log10(max(df[df$Category != "NOT", "PValue"]))
  #minq <- -log10(max(df[df$Category != "NOT", "FDR"]))
  #yhline <- min(c(minp, minq))
  if (!is.null(ymax)) {
    cat(paste0("注意：-Log2 P Value > ", ymax, "的数据被剔除!\n"))
    df <- df[-log10(df$PValue) < ymax, ]
  }
  # valcano
  p = ggplot(df, aes(x = log2FC, y = -log10(PValue))) + 
    geom_point(data = df, aes(x = log2FC, 
                              y = -log10(PValue), color = Category), 
               shape = shape, alpha = alpha)
  # mark valcano
  if (!is.null(label_geneset)) {
    df_label <- df[which(df$Gene %in% label_geneset), ]
    p <- p + 
      geom_point(data = df_label, aes(x = log2FC, y = -log10(PValue)), 
                 color = "black", size = 4) + 
      geom_point(data = df_label, aes(x = log2FC, y = -log10(PValue)), 
                 color = "white", size = 2.5) + 
      geom_point(data = df_label, aes(x = log2FC, y = -log10(PValue), 
                                      color = Category), 
                 size = 1.5) + 
      geom_text_repel(data = df_label, aes(x = log2FC, y = -log10(PValue), 
                                           label = Gene), 
                      fontface = "bold", max.overlaps = 30)
  }
  # add theme
  p <- p + theme_bw() + 
    geom_vline(xintercept = c(logfc.n, logfc.p), lty = 2) + 
    #geom_hline(yintercept = yhline, lty = 2) + 
    geom_hline(yintercept = -log10(pvalue), lty = 2) + 
    scale_x_continuous(limits = c(-xlim, xlim)) + 
    coord_fixed(ratio = (2 * xlim)/(max(-log10(df$PValue[df$PValue != 0]), 
                                        na.rm = T))) + 
    theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(color = "black")) + 
    xlab("Log2(Fold Change)") + ylab("-Log10(P Value)") + 
    scale_color_manual(values = c(NOT = color.not, 
                                  UP = color.up, DOWN = color.down))
  p
  outname <- paste0(dirclean(outdir), "/", filename.base, ".pdf")
  ggsave(outname, width = pic.w, height = pic.h)
  return(list(pic = p, data = df))
}


#' Heatmap
#' @description plot heatmap with selected gene
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
#' @param f_z number, font size
#' @param f_z.col number, col font size
#' @param angle character, col angle
#' @param cluster_rows logical, default is T
#' @param cluster_cols logical, default is F
#' @return #
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom ggsci pal_d3
#' @import RColorBrewer
#'
#' @examples #
DEGp_Pheat <- function(df, gene, pheno, outdir, f_mark='', nlevel=8,
                      pic_w = 6, pic_h = 8, c_w=1, c_h = 10, f_z = 12,
                      f_z.col= 8, angle="90", cluster_rows= T,
                      cluster_cols = F) {
  ## 得到热图数据所需格式
  #df = heat_matrix ; gene = sig_gene_100
  heat_matrix <- df[rownames(df) %in% gene,]
  # 设置分组顏色
  matdf <- pheno
  if (nlevel != 2) {
    mat_colors <- list(ggsci::pal_d3()(nlevel))
  } else {
    mat_colors <- list(c("dimgrey","darkorange"))
  }
  names(mat_colors)[1] <- names(matdf)[1]
  names(mat_colors[[1]]) <- levels(matdf[, 1])

  col <- colorRampPalette(brewer.pal(12,'Paired'))(24)
  #plot(1:24,rep(1,24),col= col,pch=16,cex=2)
  
  mkdir(dirclean(outdir))
  fname <- paste0(dirclean(outdir), '/heatmap.', f_mark, '.pdf')
  pdf(fname, width = pic_w, height = pic_h)
  #p <- pheatmap(heat_matrix, scale = 'row',cluster_rows = T,
  #              border_color = NA, show_colnames = T,
  #              show_rownames = T, fontsize = 11, fontsize_row = 7,
  #              annotation_col = matdf, annotation_colors = mat_colors,
  #              angle_col = "45",
  #              drop_levels = T, color = colp, name = 'Color',
  #              annotation_names_col = F)
  p <- pheatmap(heat_matrix, scale = "row",
                cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = colorRampPalette(colors = c("blue", "white", "red"))(100),
                annotation_col = matdf,  annotation_colors = mat_colors,
                show_colnames = T, show_rownames = T,
                border_color = NA,
                cellwidth = c_w, cellheight = c_h,
                fontsize = f_z, fontsize_col = f_z.col, 
                angle_col = angle, drop_levels = T)
  print(p)
  dev.off()
  #
  png(filename = paste0(outdir, '/heatmap.', f_mark, '.png'),
      height = pic_h, width = pic_w, units = "in", res = 1000)
  print(p)
  dev.off()
}


#' Venn Plot
#' @description Plot Venn 2-5 group data
#'
#' @param datalist list, list(diffan$deg$symbol, diffagek)
#' @param transparency number, transparency
#'
#' @return #
#' @export
#' @importFrom ggsci pal_jco
#' @importFrom scales alpha
#' @importFrom VennDiagram venn.diagram
#'
#' @examples #
DEGp_venn <- function(datalist, transparency = 0.6) {
  # datalist为list,见例子 。 注意本函数返回结果必须使用grid.draw()函数打印
  # transparency 为 0-1，表示透明度，或设置为NULL，表示不至于区域颜色
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
