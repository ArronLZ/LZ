#' @title DEGp_prepareVolcano
#' @description This function is used to prepare the data for the Volcano plot.
#'
#' @param resdf data.frame, the result of DEG analysis, must contain Gene, log2FC, PValue, FDR columns.
#' @param filterc character, "pvalue" or "fdr" or "both", which column to use as filter.
#' @param pvalue numeric, the pvalue threshold.
#' @param fdr numeric, the fdr threshold.
#' @param logfc.p numeric, the positive log2FC threshold.
#' @param logfc.n numeric, the negative log2FC threshold.
#' @param ymax numeric, the y axis threshold, the data greater ymax will be delete.
#'
#' @importFrom dplyr select
#' @return data.frame, the data for the Volcano plot.
#' @author Jiang
DEGp_prepareVolcano <- function(resdf, filterc = "pvalue", 
                                pvalue=0.05, fdr=0.1, 
                                logfc.p=1, logfc.n=-1, ymax = NULL) {
  if (!all(c("Gene", "log2FC", "PValue", "FDR") %in% names(resdf) )) {
    stop("数据表格必须包含Gene, log2FC, PValue, FDR这四列")
  }
  df <- resdf %>% dplyr::select(Gene, log2FC, PValue, FDR)
  
  # add category
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
  
  # Delete the data that is greater than ymax
  if (!is.null(ymax)) {
    cat(paste0("注意：-Log2 P Value > ", ymax, "的数据被剔除!\n"))
    df <- df[-log10(df$PValue) < ymax, ]
  }
  return(df)
}


#' DEGp_Volcano_cor
#' @description This function is used to plot the Volcano plot.
#' 
#' @param resdf data.frame, the result of DEG analysis, must contain Gene, log2FC, PValue, FDR columns.
#' @param filterc character, "pvalue" or "fdr", which column to use as filter.
#' @param pvalue numeric, the pvalue threshold.
#' @param fdr numeric, the fdr threshold.
#' @param logfc.p numeric, the positive log2FC threshold.
#' @param logfc.n numeric, the negative log2FC threshold.
#' @param ymax numeric, the y axis threshold, the data greater ymax will be delete.
#' @param ylimup numeric, the y axis limit up.
#' @param label_geneset character(n), the gene set to label.
#' @param label_move numeric, the move distance of the label.
#' @param outdir character, the output directory.
#' @param filename.base character, the output file name.
#' @param xlab_text character, the x axis label.
#' @param pic.w numeric, the width of the plot.
#' @param pic.h numeric, the height of the plot.
#' @param color.not character, the color of the data that is not significant.
#' @param color.up character, the color of the data that is up-regulated.
#' @param color.down character, the color of the data that is down-regulated.
#' @param size numeric, the size of the point.
#' @param shape numeric, the shape of the point.
#' @param alpha numeric, the transparency of the point.
#'
#' @import ggplot2
#' @return list, the list contains the plot and the data.
#' @export
#'
#' @author Jiang
DEGp_Volcano_cor <- function(resdf, filterc = "pvalue", 
                             pvalue=0.05, fdr=0.1, logfc.p=1, logfc.n=-1, 
                             ymax = NULL, ylimup = 5,
                             label_geneset = NULL, label_move = 0.2,
                             outdir="result", filename.base = "DEG_xx", 
                             xlab_text = "Pearson correlation(Pearson test)",
                             pic.w = 7, pic.h =7,
                             color.not = "#BEBEBE", color.up = "#AC001D", 
                             color.down = "#3D71AD", size=1.5, 
                             shape=16, alpha = 0.9) {
  df <- DEGp_prepareVolcano(resdf, filterc, pvalue, fdr, logfc.p, logfc.n, ymax)
  xlim = max(abs(df$log2FC))
  ylimp = max(-log10(df$PValue[df$PValue != 0]), na.rm = T) + ylimup
  ylimq = max(-log10(df$FDR[df$FDR != 0]), na.rm = T) + ylimup
  
  # plot main() ----------------------------
  # 1. plot mode(which column to use as y axis)
  if (filterc == "pvalue") {
    # y axis is -log10(PValue)
    p = ggplot(df, aes(x = log2FC, y = -log10(PValue))) + 
      geom_point(data = df, aes(x = log2FC, 
                                y = -log10(PValue), color = Category), 
                 shape = shape, alpha = alpha) +
      geom_hline(yintercept = -log10(pvalue), lty = 2) +
      coord_fixed(ratio = (2 * xlim)/(max(-log10(df$PValue[df$PValue != 0]), 
                                          na.rm = T))) +
      scale_y_continuous(limits = c(0, ylimp)) +
      ylab("-Log10(P Value)")
  } else {
    # y axis is -log10(FDR)
    p = ggplot(df, aes(x = log2FC, y = -log10(FDR))) + 
      geom_point(data = df, aes(x = log2FC, 
                                y = -log10(FDR), color = Category), 
                 shape = shape, alpha = alpha) +
      geom_hline(yintercept = -log10(fdr), lty = 2) +
      coord_fixed(ratio = (2 * xlim)/(max(-log10(df$FDR[df$FDR != 0]), 
                                          na.rm = T))) +
      scale_y_continuous(limits = c(0, ylimq)) +
      ylab("-Log10(FDR)")
  }
  
  # 2. add label
  if (!is.null(label_geneset)) {
    df_label <- df[which(df$Gene %in% label_geneset), ]
    df_label$log2FC <- ifelse(df_label$log2FC < 0, 
                              df_label$log2FC - label_move, df_label$log2FC + label_move)
    if (filterc == "pvalue") {
      p <- p + geom_text(data = df_label, 
                         aes(x = log2FC, y = -log10(PValue), label = Gene),
                         size = 2)
    } else {
      p <- p + geom_text(data = df_label, 
                         aes(x = log2FC, y = -log10(FDR), label = Gene),
                         size = 2)
    }
  }
  
  # 3. add theme
  p <- p + theme_bw() + 
    geom_vline(xintercept = c(logfc.n, logfc.p), lty = 2) + 
    scale_x_continuous(limits = c(-xlim-0.3, xlim+0.3)) + 
    theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(color = "black")) + 
    xlab(xlab_text) + 
    scale_color_manual(values = c(NOT = color.not, 
                                  UP = color.up, DOWN = color.down))
  
  # save plot
  mkdir(dirclean(outdir))
  outname <- paste0(dirclean(outdir), "/", filename.base, ".pdf")
  ggsave(outname, width = pic.w, height = pic.h)
  return(list(pic = p, data = df))
}


#' Valcano2
#' @description This function is used to plot the Volcano plot.
#' 
#' @param resdf data.frame, the result of DEG analysis, must contain Gene, log2FC, PValue, FDR columns.
#' @param filterc character, "pvalue" or "fdr", which column to use as filter.
#' @param pvalue numeric, the pvalue threshold.
#' @param fdr numeric, the fdr threshold.
#' @param logfc.p numeric, the positive log2FC threshold.
#' @param logfc.n numeric, the negative log2FC threshold.
#' @param label_geneset character(n), the gene set to label.
#' @param ymax numeric, the y axis threshold, the data greater ymax will be delete.
#' @param outdir character, the output directory.
#' @param filename.base character, the output file name.
#' @param xlab_text character, the x axis label.
#' @param pic.w numeric, the width of the plot.
#' @param pic.h numeric, the height of the plot.
#' @param color.not character, the color of the data that is not significant.
#' @param color.up character, the color of the data that is up-regulated.
#' @param color.down character, the color of the data that is down-regulated.
#' @param size numeric, the size of the point.
#' @param shape numeric, the shape of the point.
#' @param alpha numeric, the transparency of the point.
#'
#' @return list, the list contains the plot and the data.
#' @export
#' @import ggrepel
#' @include tool.R
#' @author Jiang
DEGp_Volcano2 <- function(resdf, filterc = "pvalue", 
                          pvalue=0.05, fdr=0.1, logfc.p=1, logfc.n=-1, 
                          label_geneset = NULL, 
                          ymax = NULL, 
                          outdir="result", filename.base = "DEG_xx", 
                          xlab_text = "log2(Fold Change)",
                          pic.w = 7, pic.h =7,
                          color.not = "#BEBEBE", color.up = "#AC001D", 
                          color.down = "#3D71AD",  size=1.5, 
                          shape=16, alpha = 0.9) {
  if (max(resdf$log2FC) > 30) {
    cat("警告，严重警告：基因的变化倍数列数据最大值大于30,意味着该数据可能尚未log2化\n")
  }
  df <- DEGp_prepareVolcano(resdf, filterc, pvalue, fdr, logfc.p, logfc.n, ymax)
  xlim = max(abs(df$log2FC))
  ylimp = max(-log10(df$PValue[df$PValue != 0]), na.rm = T)
  ylimq = max(-log10(df$FDR[df$FDR != 0]), na.rm = T)
  
  # plot main() ----------------------------
  # 1. plot mode(which column to use as y axis)
  if (filterc == "pvalue") {
    # y axis is -log10(PValue)
    p = ggplot(df, aes(x = log2FC, y = -log10(PValue))) + 
      geom_point(data = df, aes(x = log2FC, 
                                y = -log10(PValue), color = Category), 
                 shape = shape, alpha = alpha) +
      geom_hline(yintercept = -log10(pvalue), lty = 2) +
      coord_fixed(ratio = (2 * xlim)/(max(-log10(df$PValue[df$PValue != 0]), 
                                          na.rm = T))) +
      scale_y_continuous(limits = c(0, ylimp)) +
      ylab("-Log10(P Value)")
  } else {
    # y axis is -log10(FDR)
    p = ggplot(df, aes(x = log2FC, y = -log10(FDR))) + 
      geom_point(data = df, aes(x = log2FC, 
                                y = -log10(FDR), color = Category), 
                 shape = shape, alpha = alpha) +
      geom_hline(yintercept = -log10(fdr), lty = 2) +
      coord_fixed(ratio = (2 * xlim)/(max(-log10(df$FDR[df$FDR != 0]), 
                                          na.rm = T))) +
      scale_y_continuous(limits = c(0, ylimq)) +
      ylab("-Log10(FDR)")
  }
  
  # 2. add label
  if (!is.null(label_geneset)) {
    df_label <- df[which(df$Gene %in% label_geneset), ]
    if (filterc == "pvalue") {
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
                        fontface = "bold", max.overlaps = 3000)
    } else {
      p <- p + 
        geom_point(data = df_label, aes(x = log2FC, y = -log10(FDR)), 
                   color = "black", size = 4) + 
        geom_point(data = df_label, aes(x = log2FC, y = -log10(FDR)), 
                   color = "white", size = 2.5) + 
        geom_point(data = df_label, aes(x = log2FC, y = -log10(FDR), 
                                        color = Category), 
                   size = 1.5) + 
        geom_text_repel(data = df_label, aes(x = log2FC, y = -log10(FDR), 
                                             label = Gene), 
                        fontface = "bold", max.overlaps = 3000)
    }
  }
  
  # add theme
  p <- p + theme_bw() + 
    geom_vline(xintercept = c(logfc.n, logfc.p), lty = 2) + 
    scale_x_continuous(limits = c(-xlim-0.3, xlim+0.3)) + 
    theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(color = "black")) + 
    xlab(xlab_text) + 
    scale_color_manual(values = c(NOT = color.not, 
                                  UP = color.up, DOWN = color.down))
  
  mkdir(dirclean(outdir))
  outname <- paste0(dirclean(outdir), "/", filename.base, ".pdf")
  ggsave(outname, width = pic.w, height = pic.h)
  return(list(pic = p, data = df))
}
