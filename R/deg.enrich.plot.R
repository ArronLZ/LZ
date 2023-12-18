#' Prepare Dotplot
#' @description Prepare Dotplot
#' 
#' @param df dataframe, from go or kegg analysis
#' @param head numeric or NULL, head = 20：filter the head 20 row of df, head = NULL：all of df 
#' @param delete numeric or NULL, delete = c(4, 7)：delete 4, 7, NULL: not delete
#'
#' @return #
#' @export
#'
#' @examples # 
#' @author Jiang
DEGp_prepareDotplot <- function(df, head = 20, delete = NULL) {
  # df为dataframe, 为go或kegg分析后取出的表格数据
  #   格式如下（需有Description, pvalue, Count, GeneRatio列）
  # 建议事先将df排序好或者挑选好
  # head = 20, 表示取表格前20, head = NULL, 表示取全部 
  # delete = c(4, 7)， 表示删除表格中的第4和7行，delete = NULL, 表示不删除表格数据
  df$GeneRatio <- strsplit(df$GeneRatio, '/') %>% 
    lapply(., function(x) round(as.numeric(x[1])/as.numeric(x[2]), 2)) %>% 
    as.numeric()
  if (!is.null(head)) {
    df <- df %>% head(head) %>% na.omit()
  }
  if (!is.null(delete)) {
    df <- df[-delete, ]
  }
  return(df)
}


#' Dotplot 
#' @description Dotplot of go or kegg result(DEGp_prepareDotplot result)
#' 
#' @param df dataframe, from DEGp_prepareDotplot result
#' @param title character
#' @param resultdir cahracter, picture save dir
#' @param filemark character, picture filename mark
#' @param pic.save logical, save picutre or not?
#' @param pw number, picuter width
#' @param ph number, picuter height
#'
#' @return # 
#' @export
#'
#' @examples #
#' @author Jiang
DEGp_Dotplot <- function(df, title='xxx', resultdir, filemark, pic.save =T, 
                         pw=8, ph=10) {
  # data格式如下（需有Description,pvalue,Count,GeneRatio列）
  #                             Description       pvalue Count GeneRatio
  #GO:0023061                signal release 5.246744e-06    29      0.06
  #GO:0050673 epithelial cell proliferation 6.129368e-06    27      0.06
  #GO:0045444      fat cell differentiation 7.622508e-06    18      0.04
  # resultdir = './result_stringtie/p005fc15'
  # filemark = 'GO_BP_top20_common'
  # 此函数可能需要使用scale_x_continuous调整x轴刻度
  dotplot <- ggplot(cbind(df, Order = nrow(df):1)) +
    geom_point(mapping = aes(x = -log10(pvalue), y = Order, 
                             size = GeneRatio, fill = qvalue),
               shape = 21) + 
    scale_fill_gradientn(colours = c("grey", "gold", "red")) + #自定义配色
    scale_y_continuous(position = "left", 
                       breaks = 1:nrow(df), 
                       labels = Hmisc::capitalize(rev(df$Description))) +
    #scale_x_continuous(breaks = c(4,5,6,7),
    #                   #breaks = seq(0, xmax+5, 5),
    #                   limits = c(4,7),
    #                   expand = expansion(mult = c(.05, .05))) + #两边留空
    labs(x = "-log10(P-value)", y = NULL) +
    guides(size = guide_legend(title = "GeneRatio"),
           fill = guide_colorbar(title = "qvalue")) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtheme.update.text()
  if (pic.save == T) {
    dotplot %>% ggplotGrob() %>% cowplot::plot_grid()
    fname=paste0(resultdir, '/', filemark,'.pdf')
    ggsave(fname, width = pw, height = ph)
  }
  return(dotplot)
}


#' Dotplot2
#' @description Dotplot of go or kegg result(clusterProfiler enrich result) or your own data\cr
#' the dataframe must contain colname : "Description", "GeneRatio", "pvalue", "qvalue", "Count"\cr
#' note: the col GeneRatio data type should be numeric, if not please trans to numeric \cr
#' this function do not need prepare data  using DEGp_prepareDotplot() 
#' @param df data.frame the dataframe must contain colname : "Description", "GeneRatio", "pvalue", "qvalue", "Count"
#' @param head numeric deafault is 20：keep the head 20 row of df, if not set, use all of the df
#' @param delete numeric or NULL delete = c(4, 7)：delete 4, 7, NULL: not delete
#' @param method number picture plot method, can be set one of 1:4, default is 1
#' @param title character the picture title
#' @param resultdir character outdir name where the picture save
#' @param filemark character the picture filename, default is "xxx"
#' @param pic.save logical save the picture or not,, default is T
#' @param pw number the picture width, default is 8
#' @param ph number the picture height, default is 10
#' @param pcompress number create a blank area at the top of the picture, default is 1
#' @param size1 number the dot size start, default is 3
#' @param size2 number the dot size end, default is 5
#' @param color1 character color rgb
#' @param color2 character color rgb
#'
#' @return ggplot2 picutre
#' @export
#'
#' @author Jiang
DEGp_Dotplot2 <- function(df, head = 20, delete = NULL,
                          method = 1, title = "xxx", 
                          resultdir, filemark, pic.save = T, 
                          pw = 8, ph = 10, pcompress = 1, 
                          size1=3, size2=5, 
                          color1="#AC001D", color2="#888888") {
  # condition check
  stopifnot("Description" %in% names(df))
  stopifnot("GeneRatio" %in% names(df))
  stopifnot("pvalue" %in% names(df))
  stopifnot("qvalue" %in% names(df))
  stopifnot("Count" %in% names(df))
  # tryCatch({
  #   stopifnot(3 > 4)
  # }, error = function(e) {
  #   message("x必须大于y，请检查。")
  # })
  
  # data clean ------------------
  predata <- function(df, head = 20, delete = NULL) {
    if(!is.numeric(df$GeneRatio)) {
      # only in the condition that the GeneRatio if result of clusterProfiler,
      #  its GeneRatio is two number sperated by "/"
      df$GeneRatio <- strsplit(df$GeneRatio, "/") %>% 
        lapply(., function(x) round(as.numeric(x[1])/as.numeric(x[2]), 2)) %>% 
        as.numeric()
    }
    if (!missing(head)) {
      df <- df %>% head(head) %>% na.omit()
    }
    if (!is.null(delete)) {
      df <- df[-delete, ]
    }
    return(df)
  }
  df <- predata(df, head = head, delete = delete)
  
  # plot picture ---------------
  p <- ggplot(cbind(df, Order = nrow(df):1))
  if (method == 1) {
    p <-  p + 
      geom_point(mapping = aes(x = -log10(pvalue), y = Order, 
                               size = GeneRatio, fill = qvalue), shape = 21)
  } else if (method == 2) {
    p <- p + 
      geom_point(mapping = aes(x = GeneRatio, y = Order, 
                               size = pvalue, fill = qvalue), shape = 21)
  } else if (method == 3) {
    p <- p + 
      geom_point(mapping = aes(x = -log10(pvalue), y = Order, 
                               size = Count, fill = pvalue), shape = 21)
  } else if (method == 4) {
    p <- p + 
      geom_point(mapping = aes(x = -log10(qvalue), y = Order, 
                               size = Count, fill = qvalue), shape = 21)
  } else {
    p <- p + 
      geom_point(mapping = aes(x = -log10(pvalue), y = Order, 
                               size = GeneRatio, fill = qvalue), shape = 21)
  }
  p <- p + scale_y_continuous(position = "left", # expand = c(0.05,0.5),
                              breaks = 1:nrow(df),
                              labels = Hmisc::capitalize(rev(df$Description))) +
    labs(y = NULL) +
    scale_size(range = c(size1, size2)) +
    scale_fill_gradientn(colours = c(color1, color2)) + 
    ggtitle(paste0(paste(rep("\n", pcompress), collapse = ""), title)) +
    theme_classic() +
    # theme_bw() + 
    theme(panel.grid = element_blank()) +
    ggtheme.update.text() +
    guides(# size = guide_legend(reverse = TRUE),
           fill = guide_colourbar(reverse = TRUE))
  p # %>% ggplotGrob() %>% cowplot::plot_grid()
  
  # data output --------------
  if (pic.save == T) {
    fname = paste0(resultdir, "/", filemark, ".pdf")
    pdf(fname, width = pw, height = ph)
    print(p)
    dev.off()
  }
  return(p)
}


#' Prepare data for DEGp_GOALL_barplot
#' @description Prepare data for DEGp_GOALL_barplot. Take top n number of each type in the dataframe.
#' @param data dataframe, GO ALL dataframe
#' @param n_type charater, a colname of data, the type of GO ALL dataframe, such as n_type='ONTOLOGY'
#' @param n_qu character, a colname of data, sort by this param, such as n_qu='pvaule'
#' @return #
#' @export
#'
#' @examples #
#' @author Jiang
DEGp_prepareGOALL <- function(data, n_type='ONTOLOGY', n_qu='pvaule') {
  ## 自定义函数，提取每个分类的前10名
  ## data=up,为以数据框，其含有名为ncol的一列；
  ## n_type='ONTOLOGY', n_qu='pvaule'为data中的列名
  # 此函数即为，取data$ONTOLOGY的每一类pavlue列的前10名(最小的10个pvalue)
  lapply(names(data[, n_type] %>% table), function(x) 
    data[data[,n_type] == x,] %>% arrange(n_qu) %>% .[1:10,]) %>% 
    do.call(rbind, .) %>% na.omit()
}


#' GO all barplot
#' @description GO all barplot
#' @param GOlist.df dataframe, GO result
#' @param resultdir character, the picture outdir
#' @param filemark character, the picture filename mark
#' @param pic.save logic, save picture or not
#' @return #
#' @export
#'
#' @author Jiang
DEGp_GOALL_bar <- function(GOlist.df, resultdir, filemark, pic.save = T) {
  # GOlist为jGO自定义函数返回对象，该对象为clusterProfiler::enrichGO返回对象的list集合
  # GOlist.df为GOlist$ALL@result, 为一个data.frame
  # 或者可更改第一句，输入数据框含ONTOLOGY,pvalue,Description列
  # up_df <- GOlist$ALL@result
  up_df <- GOlist.df
  up <- up_df %>% DEGp_preparGOALL(., 'ONTOLOGY', 'pvalue')
  up <- up %>% arrange(ONTOLOGY, desc(pvalue))
  dy_levels <- c(filter(up, ONTOLOGY=='MF')[,"Description"],
                 filter(up, ONTOLOGY=='CC')[,"Description"],
                 filter(up, ONTOLOGY=='BP')[,"Description"])
  up$Description <- factor(up$Description, levels = dy_levels)
  id_up <- levels(up$Description)
  # GO
  #display.brewer.all(colorblindFriendly = T) #展示所有的色板
  #display.brewer.pal(12,'Paired')# 展示'Accent'色板中8个颜色
  col <- colorRampPalette(brewer.pal(12,'Paired'))(24)
  #plot(1:24,rep(1,24),col= col,pch=16,cex=2)
  colp <- colorRampPalette(c(col[10],col[6],col[3]))(100)
  #plot(1:10,col=col[6], pch=16, cex=2)
  p_up <- ggplot(up, aes(Description, -log(pvalue, 10), fill=ONTOLOGY)) +
    geom_col(color = 'black', width = 0.6) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent')) +
    theme(axis.line.x = element_line(colour = 'black'), 
          axis.line.y = element_line(colour = 'transparent'), 
          axis.ticks.y = element_line(colour = 'transparent'),
          axis.text = element_text(face = 'plain', size = 12)) +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
    coord_flip() +
    geom_hline(yintercept = 0) +
    labs(x = '', y = '', title = 'UP') +
    scale_y_continuous(expand = c(0, 0), breaks = c(-15, -10, -5, 0), 
                       labels = as.character(c(-15, -10, -5, 0))) +     #这儿更改间距设置
    scale_x_discrete(labels = id_up)
  if (pic.save == T) {
    fname=paste0(resultdir, '/', filemark,'.pdf')
    ggsave(fname, width = 6, height = 10)
  }
  return(p_up)
}
