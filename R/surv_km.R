#' @title ostime cut using different methods 
#' @description Use different methods to cut the data into two groups and then fuse this grouping information into data
#' 
#' @param eset_os.df data.frame, the data.frame of the clinical and gene expression data
#' @param tagetGene character, the gene name, the column name of the gene expression data
#' @param method character, default "mean", the method to cut the gene expression data, \cr
#' can be "mean", "upper75 vs low25", "median", "upper25 vs low75", "auto cut point", "upper25 vs low25" or a number
#' @importFrom survminer surv_cutpoint
#' 
#' @return a list
#'
#' @examples
#' \dontrun{
#' eset_os.df
#'      OS.time OS  EGFR  TP53 ...
#' P1     3450  0    4.6   6.6 ...
#' P2     1280  1    1.0   6.0 ...
#' P3      455  1    3.0   4.0 ...
#' }
Surv_ostime_cut <- function(eset_os.df, tagetGene, method="median") {
  stopifnot(is.data.frame(eset_os.df), tagetGene %in% colnames(eset_os.df))
  if (is.numeric(method)) {
    cat(method, "\n")
  } else {
    stopifnot(method %in% c("mean", "upper75 vs low25", "median", 
                            "upper25 vs low75", "auto cut point", 
                            "upper25 vs low25"))
  }
  df <- eset_os.df
  tagetGene.five <- fivenum(df[,tagetGene])
  if (method == "mean") {
    # mean
    cut.p <- mean(df[,tagetGene])
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, "h.mean","l.mean"))
  } else if (method == "upper75 vs low25") {
    # upper 75 vs low 25
    cut.p <- tagetGene.five[2]
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, "high75","low25"))
  } else if (method == "median") {
    # median
    cut.p <- tagetGene.five[3]
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, "h.median","l.median"))
  } else if (method == "upper25 vs low75") {
    # upper 25 vs low 75
    cut.p <- tagetGene.five[4]
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, "high25","low75"))
  } else if (method == "auto cut point") {
    # cut poion
    res.cut <- surv_cutpoint(df, #数据集
                             time = "OS.time", #生存状态
                             event = "OS", #生存时间
                             variables = tagetGene) #需要计算的数据列名
    cut.p <- res.cut$cutpoint$cutpoint                         
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, "h.point","l.point"))
  } else if (method == "upper25 vs low25") {
    cut.p <- c(tagetGene.five[2], tagetGene.five[4])
    df.low25 <- df[df[,tagetGene] <= cut.p[1], ]
    df.low25$Group <- "low25"
    df.high25 <- df[df[,tagetGene] > cut.p[2], ]
    df.high25$Group <- "high25"
    df <- rbind(df.low25, df.high25) %>% data.frame(check.names = F)
  } else if (is.numeric(method)) {
    cut.p <- method
    df$Group <- as.vector(ifelse(df[,tagetGene] > cut.p, 
                                 paste0("h.", cut.p), paste0("l.", cut.p)))
  }
  return(list(eset_clin_p=df, cut.p=cut.p, five=tagetGene.five))
}


#' Survival analysis and plot
#' @description Survival analysis and plot using existed group information
#'
#' @param df data.frame, the data.frame of the clinical and gene expression data
#' @param col_name character, the column name of the group column, default is "Group"
#' @param num.tran numeric, the number to transform the OS.time
#' @param break.x numeric, the break x
#' @param main.text character, the main text of the plot
#' @param color character, the color of the plot
#' @importFrom survminer ggsurvplot
#' @importFrom survival survdiff survfit
#'
#' @return a list
#' @export
#'
#' @examples
#' \dontrun{
#' df
#'      OS.time OS  EGFR  Group 
#' P1     3450  0    4.6   High   
#' P2     1280  1    1.0    Low   
#' P3      455  1    3.0    Low   
#' res <- Surv_diffG(df, col_name, num.tran=365, break.x=2.5, main.text,
#'                  color=c("#eaa8a8", "#729fc9"))
#' res <- Surv_diffG(df, col_name, num.tran=365, break.x=2.5, main.text,
#'                  color=c("#eaa8a8", "#729fc9", "#f6f1af"))
#' }
Surv_diffG <- function(df, col_name, num.tran=365, break.x=2.5, main.text,
                      color=c("#eaa8a8", "#729fc9")
                      ) {
  stopifnot(is.data.frame(df), !missing(col_name), col_name %in% colnames(df))
  # 指定分组信息为Group
  if (col_name != "Group") {
    names(df)[names(df) == col_name] = "Group"
  }
  df$OS.time <- df$OS.time/num.tran
  diff <- survdiff(Surv(OS.time, OS) ~ Group, data = df)
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = df)
  pValue <- 1 - pchisq(diff$chisq, 
                       df = (length(levels(as.factor(df$Group))) - 1)
  )
  pValue <- signif(pValue, 4)
  # plot
  #"#4DAF4A"
  if (length(levels(as.factor(df$Group))) == length(color)) {
    mycolor <- color
  } else {
    stop("The color length is not equal to the group length, use default color\n")
  }
  ggs <- ggsurvplot(fit, # 创建的拟合对象
                    data = df,  # 指定变量数据来源
                    palette = mycolor,
                    conf.int = F, # 显示置信区间
                    pval = T, # 添加P值
                    surv.median.line = "hv",  # 添加中位生存时间线
                    risk.table = TRUE, # 添加风险表
                    xlab = "Time in Years", # 指定x轴标签
                    legend = c(0.8,0.8), # 指定图例位置
                    legend.title = main.text, # 设置图例标题
                    #legend.labs = 1:5, # 指定图例分组标签
                    tables.y.text = F,
                    #break.time.by = 1,
                    break.x.by = break.x, # 设置x轴刻度间距
                    ggtheme = theme_survminer(font.legend = c(11, "plain", "black")),
                    tables.theme = theme_survminer(font.main = 11) )
  ggs$plot <- ggs$plot + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0))
  return(list(pValue=pValue, fit=fit, pic=ggs))
}


#' Survival analysis and plot
#' @description Survival analysis and plot using multi methods
#' @param eset_os data.frame, the data.frame of the clinical and gene expression data
#' @param tagetGene character, the gene name
#' @param method character, the method to cut the gene expression data, \cr
#' can be "mean", "upper75 vs low25", "median", "upper25 vs low75", "auto cut point", "upper25 vs low25" or a number
#' @param num.tran number, the number to transform the OS.time, mean: OS.time/num.tran
#' @param break.x number, the break x
#' @param col_name character, the column name of the group column, default is "Group"
#' @param color character, the color of the plot
#'
#' @return a list
#' @export
#'
#' @examples
#' \dontrun{
#' eset_os
#'      OS.time OS  EGFR  TP53 ...
#' P1     3450  0    4.6   6.6 ...
#' P2     1280  1    1.0   6.0 ...
#' P3      455  1    3.0   4.0 ...
#' res <- Surv_diffONE(eset_os, tagetGene="EGFR", method = "median",
#'                     num.tran=365, break.x=2.5, col_name="Group",
#'                     color=c("#eaa8a8", "#729fc9"))
#' }
Surv_diffONE <- function(eset_os, tagetGene, method = "median",
                           num.tran=365, break.x=2.5, col_name="Group",
                           color=c("#eaa8a8", "#729fc9")) {
  eset_os.cut <- Surv_ostime_cut(eset_os.df = eset_os, tagetGene = tagetGene, 
                                 method = method)
  ps <- Surv_diffG(eset_os.cut$eset_clin_p, main.text = tagetGene, col_name,
                  num.tran, break.x, color)
  return(list(eset_os.cut = eset_os.cut, pic = ps))
}