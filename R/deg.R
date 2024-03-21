#' RNAseq differential analysis DESeq2 001<dds>
#' @description The function is designed to perform differential analysis by DESeq2<dds>.
#' @param exprset.group the list consisted of eset, group and mark.
#' @param batch batch number
#'
#' @return dds ( DESeq(dds, parallel = T) )
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @import BiocParallel
#' 
#' @author Jiang
#' @examples
#' # Run: diffan <- DEG_DESeq2.dds(exprset.group) # batch = F
#' # exprset.group 由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
#' # ** 如果bacth=T，则group中第二列为batch列。
#' # 1. eset：
#' # ** 表达谱数据格式(数值型，整型)
#' #       | sample1 | sample2 | sample3  | sample4
#' # gene1 |  34     |  23     |  56      |  23
#' # gene2 |  35     |  23     |  12      |  23
#' # gene3 |  12     |  78     |  78      |  78
#' # 2. group：
#' # ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eseti)
#' # BATCH列可以没有(如果设置了bacth=T，BATCH列必须有，且必须位于第二列)
#' #            Type     BATCH
#' # sample1  | tumor   | 1
#' # sample2  | tumor   | 1
#' # sample3  | control | 2
#' # sample4  | control | 2
#' # 3. f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
#' # 默认的比法是后比前（字符排序在后面的比前面的，如示例中的就是 t vs c,因为c排在前面）
#' 
#' # LZ包将DESeq2差异分析整合拆解成了三个函数，因此需要如下代码才能完整进行差异分析。
#' # 差异分析 deseq2三部曲
#' # dds <- DEG_DESeq2.dds(exprset.group=glist, batch = F)  # 差异分析1：准备工作
#' # DEG_DESeq2.pca(dds, outdir = outdirsub)   # 画PCA图。 此处有warning信息，不用管。
#' # dds_list <- DEG_DESeq2.ana(dds)           # 差异分析1：得到结果并进行初步筛选
DEG_DESeq2.dds <- function(exprset.group, batch = F) {
  #keepGene <- rowSums(edgeR::cpm(exprset)>0) >=2
  #exprset <- exprset[keepGene,]
  # 1.1 group
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  if (batch == T) {
    design_type <- as.formula(paste0('~ ', colnames(pheno)[1],
                                     '+', colnames(pheno)[2]))
  }
  # 不考虑批次效应才需要下列两行代码
  if (batch == F) {
    pheno <- data.frame(Type=pheno[,1])
    design_type <- as.formula(paste0('~ ', colnames(pheno)[1]))
  }
  # 1.2 dds
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = exprset,
                                colData = pheno,
                                design = design_type)
  keep <- rowSums(counts(dds) >= 10) >= ncol(exprset)/3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数
  dds <- dds[keep, ]
  # 1.3
  dds <- DESeq2::DESeq(dds, parallel = T)
  # ----------------------------------
  # DEG_DESeq2.pca(deseq2.dds=dds, outdir = outdir)
  # DEG_DESeq2.ana
  return(dds)
}


#' RNAseq differential analysis DESeq2 002<pca>
#' @description The function is designed to perform differential analysis by DESeq2<pca>.
#' @param dds the result of DEG_DESeq2.dds
#' @param outdir output dir, pca plot output dir, if the outdir exit the function will remove it
#' @return no return, pca plot to disk
#'
#' @export
#' @importFrom DESeq2 vst
#' @importFrom DESeq2 plotPCA
#' @importFrom DESeq2 plotDispEsts
#' @importFrom factoextra hcut
#' @importFrom factoextra fviz_dend
#'
#' @author Jiang
DEG_DESeq2.pca <- function(dds, outdir) {
  # if (!dir.exists(outdir)) {
  #   dir.create(outdir, recursive = T, showWarnings = T)
  # }
  outdir <- dirclean(outdir)
  mkdir.p(outdir)
  ## 数据转换
  # library(factoextra)
  vsd <- DESeq2::vst(dds, blind=FALSE)
  plotpca <- DESeq2::plotPCA(vsd, intgroup="Type")
  sampleDists <- stats::dist(t(assay(vsd)))
  h.data <- factoextra::hcut(sampleDists, k = 2, stand = FALSE, 
                             hc_method = "average")
  ### plot
  pdf(paste0(outdir, '/pre_PCA&Dispersion.plot', '.pdf'), width = 8, height = 7)
  # 1. Dispersion plot
  DESeq2::plotDispEsts(dds, main="Dispersion Plot")
  pf <- factoextra::fviz_dend(h.data, rect_fill = T, cex = 0.8, color_labels_by_k = T, horiz = T)
  print(pf)
  # 2. PCA plot
  p <- ggplot(plotpca$data, aes(x=PC1, y=PC2, color=Type)) +
    geom_point(size=3) +
    geom_text(aes(label=rownames(plotpca$data), y=PC2+0.1), hjust=-0.1)+
    labs(title = 'PCA analysis of samples',
         x=plotpca$labels$x, y=plotpca$labels$y) + # 坐标轴标签
    scale_x_continuous(limits = range(plotpca$data$PC1) * 1.8) + # 坐标轴范围
    theme_bw() + #coord_flip() + #主题
    theme(panel.grid.minor = element_blank(), #网格及边框
          panel.grid.major = element_line(colour = 'gray', linetype = 'dashed', linewidth = 0.2),
          panel.border = element_rect(fill=NA, color = 'black', linewidth = 1)) +
    theme(axis.text = element_text(face = 'plain', size = 12),  # 坐标轴刻度字体
          axis.title = element_text(face = 'plain', size = 12), # 坐标轴标签字体
          plot.title = element_text(hjust = 0.5, face = 'bold'))
  print(p)
  dev.off()
  ##
  png(filename = paste0(outdir,  "/pre_disp", ".png"),
      height = 7, width = 8, units = "in", res = 1000)
  DESeq2::plotDispEsts(dds, main="Dispersion Plot")
  dev.off()
  png(filename = paste0(outdir,  "/pre_cluster_tree", ".png"),
      height = 7, width = 8, units = "in", res = 1000)
  print(pf)
  dev.off()
  png(filename = paste0(outdir,  "/pre_pca", ".png"),
      height = 7, width = 8, units = "in", res = 1000)
  print(p)
  dev.off()
}


#' RNAseq differential analysis DESeq2 003<dds.filter_result>
#' @description The function is designed to perform differential analysis by DESeq2<dds.filter_result>.
#' @param dds the result of DEG_DESeq2.dds
#' @param pval PValue
#' @param fdr FDR
#' @param logfc log2(fold change)
#' @return list(list(groupdata = dds@colData, resdf=res_df, deg=etSig))
#'
#' @export
#' @importFrom DESeq2 results
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' 
#' @author Jiang
DEG_DESeq2.ana <- function(dds, pval=0.05, fdr=0.1, logfc=1) {
  # resultsNames(dds)                # 查看结果的名称。
  # dds$type                         # 默认后比前
  type1 <- levels(dds$Type)[1]
  type2 <- levels(dds$Type)[2]
  # 如不指定contrast，默认第二因子比第一因子。建议将对照的因子水平放在第一个。
  # 也可以通过results(dds, contrast=c("type", type2, type1)): type2 比 type 1
  res <- DESeq2::results(dds)
  # summary(res)      #看一下结果的概要信息，p值默认小于0.1。

  ### 差异分析总表
  res_df <- data.frame(res)
  eset_norma <- data.frame(counts(dds, normalized=TRUE), check.names = F)
  eset_norma <- log2(eset_norma + 1)
  res_df <- merge(res_df, eset_norma, by="row.names", sort=FALSE)
  res_df <- res_df %>% dplyr::rename(Gene = 'Row.names',
                                     log2FC = log2FoldChange,
                                     PValue = pvalue,
                                     FDR = padj) %>%
    dplyr::arrange(desc(log2FC), PValue) %>%
    dplyr::select(Gene, log2FC, PValue, FDR, everything())
  # 差异基因筛选
  etSig <- res_df[which(
    res_df$PValue < pval & res_df$FDR < fdr & abs(res_df$log2FC) > logfc), ]
  return(list(groupdata = dds@colData, resdf=res_df, deg=etSig))
}


#' RNAseq differential analysis edgeR
#' @description The function is designed to perform differential analysis by edgeR.
#' @details Following data is requried: express, group.
#'
#' @param exprset.group x is a 3-object list: list(eset=exprset, group=phen, f_mark="")
#' @param coef default is 2, meanig the 2nd fator vs 1st factor of group.
#' @param pval default is 0.05
#' @param fdr default is 0.1
#' @param logfc default is 1
#' @return a list: list(mid = qlf, groupdata = qlf$samples, resdf = et.norm, deg = etSig)
#' qlf: mid data of edgeR, et.norm: all result of DEG, sig result of DEG.
#'
#' @export
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr rename select arrange
#' @import edgeR
#' 
#' @examples 
#' # Run: diffan <- DEG_edgeR(exprset.group) # pval=0.05, fdr=0.1, logfc=1
#' # exprset.group 由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
#' # 1. eset：
#' # ** 表达谱数据格式(数值型，整型)
#' #       | sample1 | sample2 | sample3  | sample4
#' # gene1 |  34     |  23     |  56      |  23
#' # gene2 |  35     |  23     |  12      |  23
#' # gene3 |  12     |  78     |  78      |  78
#' # 2. group：
#' # ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eseti)
#' # BATCH列可以没有
#' #            Type     BATCH
#' # sample1  | tumor   | 1
#' # sample2  | tumor   | 1
#' # sample3  | control | 2
#' # sample4  | control | 2
#' # 3. f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
#' # 默认的比法是后比前（字符排序在后面的比前面的，如示例中的就是 t vs c,因为c排在前面）
DEG_edgeR <- function(exprset.group, pval=0.05, fdr=0.1, logfc=1, coef = 2) {
  # 此函数为edegR差异分析，返回两个表格组成的list：
  #   1.全部的分析表格 + 2.差异表格(P<0.05&FDR<0.1&abs(et$log2FC)>1)
  # exprset.group要求为list (主要为了方便lapply函数做批量差异分析)

  # 1.0 count，group
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  exprset <- exprset[, rownames(pheno)]

  # 1.1 构建DEGList
  y <- edgeR::DGEList(counts = exprset, group = pheno[,1])
  keep <- edgeR::filterByExpr(y)
  #keep <- filterByExpr(y, min.count = 10, min.total.count = 15,
  #                     large.n = 10, min.prop = 0.7)
  y <- y[keep, , keep.lib.sizes=FALSE]
  # 进行TMM标准化
  y.norm <- edgeR::calcNormFactors(y, method = "TMM")
  # 获取cpm TMM值
  cpm.tmm <- edgeR::cpm(y.norm, log = F)
  cpm.tmm <- log2(cpm.tmm + 1)

  # 1.2 计算离散因子
  y <- edgeR::calcNormFactors(y)
  # design.form <- as.formula(paste("~", "group"))
  # # y$samples$group <- relevel(y$samples$group, ref="C")
  design <- stats::model.matrix(~group, data =  y$samples)
  y <- edgeR::estimateDisp(y, design)
  #To perform quasi-likelihood F-tests:
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = coef)

  # 获取结果
  # 获取排名靠前的基因，这里设置n=80000是为了输出所以基因
  et <- topTags(qlf, n = 80000)
  et <- as.data.frame(et) # 转换为数据框类型
  et$placeholder <- "Farmer"
  et.norm <- merge(et, cpm.tmm, by = 0)
  et.norm <- tibble::column_to_rownames(et.norm, var = "Row.names")
  et.norm$Gene <- rownames(et.norm)
  et.norm <- et.norm %>% dplyr::rename(log2FC=logFC) %>%
    dplyr::arrange(desc(log2FC), PValue) %>%
    dplyr::select(Gene, log2FC, PValue, FDR, everything())
  # 差异基因筛选
  etSig <- et.norm[which(et.norm$PValue < pval &
                           et.norm$FDR < fdr & abs(et.norm$log2FC) > logfc),]
  return(list(mid = qlf, groupdata = qlf$samples, resdf = et.norm, deg = etSig))
}


#' RNAseq differential analysis Limma voom
#' @description The function is designed to perform differential analysis by Limma voom.
#'
#' @param exprset.group the list consisted of eset, group and mark.
#' @param pval default 0.05
#' @param fdr default 0.1
#' @param logfc default 1
#'
#' @return list(groupdata = pheno, resdf=tab, deg=deg)
#' @export
#' @import limma
#' @importFrom  tibble column_to_rownames
#'
#' @examples
#' # # Run: diffan <- DEG_voom(exprset.group) # pval=0.05, fdr=0.1, logfc=1
#' # exprset.group 由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
#' # 1. eset：
#' # ** 表达谱数据格式(数值型，整型)
#' #       | sample1 | sample2 | sample3  | sample4
#' # gene1 |  34     |  23     |  56      |  23
#' # gene2 |  35     |  23     |  12      |  23
#' # gene3 |  12     |  78     |  78      |  78
#' # 2. group：
#' # ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eseti)
#' # BATCH列可以没有
#' #            Type     BATCH
#' # sample1  | tumor   | 1
#' # sample2  | tumor   | 1
#' # sample3  | control | 2
#' # sample4  | control | 2
#' # 3. f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
#' # 默认的比法是后比前（字符排序在后面的比前面的，如示例中的就是 t vs c,因为c排在前面）
DEG_voom <- function(exprset.group, pval=0.05, fdr=0.1, logfc=1) {
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  if ( !all(rownames(pheno) == colnames(exprset)) ) {
    stop("colnames of exprset do not equal to rowname of pheno")
  }
  
  #pheno$type
  design <- model.matrix(~0 + pheno[,1])
  colnames(design) <- levels(pheno[,1])
  colnames(design) <- c('con','trt')  # 第二因子比第一因子
  #contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
  contrast.matrix <- makeContrasts(trt-con, levels = design)  # 第二因子比第一因子
  
  # edgeR
  y <- DGEList(counts=exprset)
  keep <- filterByExpr(y)
  #keep <- filterByExpr(y, design, min.count = 15, min.total.count = 20,
  #                     large.n = 10, min.prop = 1)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # 进行TMM标准化 ------
  y.norm <- calcNormFactors(y, method = "TMM")
  # 获取cpm TMM值 
  cpm.tmm <- cpm(y.norm, log = F)
  cpm.tmm <- log2(cpm.tmm + 1)
  # / ---------
  
  y <- calcNormFactors(y)
  # voom
  v <- voom(y, design, plot=F)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  #topTable(fit, coef=ncol(design))
  tab <- topTable(fit, sort.by = "P", number = Inf)
  tab <- cbind(rownames(tab), tab)  # 将行名粘贴为数据框的第一列
  colnames(tab) <- c("Gene", "log2FC", "AveExp", "T", "PValue", "FDR", "B")
  tab <- merge(tab, cpm.tmm, by = 0)
  tab <- tibble::column_to_rownames(tab, var = "Row.names")
  # tab <- tab %>% dplyr::select(Gene, log2FC, PValue, FDR, everything())
  tab <- tab %>% dplyr::arrange(desc(log2FC), PValue)
  
  deg <- tab[which(tab$PValue < pval & tab$FDR < fdr & abs(tab$log2FC) > logfc),]
  return(list(groupdata = pheno, resdf=tab, deg=deg))
}


#' General compare the differnt ID(Gene) between the group B vs A
#' @description General compare the differnt ID(Gene) between the group B vs A
#'
#' @param eset data.frame, the exprs data, must be log2 transformed
#' @param group data.frame, the pheno information, the first coloumn is the type information and must be factor.
#' @param pval number, the pvalue cutoff
#' @param fdr number, the fdr cutoff
#' @param logfc number, the log2(fc) value cutoff
#'
#' @return # all deg dataframe(resdf)
#' @export
#' @importFrom limma lmFit eBayes topTable
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' pn_eset <- data.table::fread("./protein_count_pn.csv", data.table = F)
#' pn_eset <- pn_eset %>% tibble::column_to_rownames(var = "Protein")
#' 
#' group <- read.csv("./group.csv", row.names = 1)
#' group$Type <- as.factor(group$Type)
#' 
#' df <- log2(pn_eset + 1)
#' diff <- limma.general(eset = df, group = group)
#' }
limma.general <- function(eset, group, pval = 0.05, fdr = 0.1, logfc = log2(2)) {
  names(group)[1] <- "Type"
  stopifnot(max(eset) < 25, all(colnames(eset) == rownames(group)),
            is.factor(group$Type))
  group_list <- group$Type
  # ~法,截距法
  design <- model.matrix(~group_list) 
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(eset)
  # ~三部曲
  fit <- lmFit(eset, design)
  fit <- eBayes(fit, trend = TRUE)
  tab <- topTable(fit, coef = 2, number = Inf)
  tab <- merge(tab, eset, by = 0)
  tab <- tab %>%
    dplyr::rename(Gene=Row.names, log2FC=logFC, PValue=P.Value, FDR=adj.P.Val) %>%
    # dplyr::select(Gene, log2FC, PValue, FDR, names(eset), everything()) %>%
    dplyr::arrange(desc(log2FC), PValue)
  
  deg <- tab[which(tab$PValue < pval & tab$FDR < fdr & abs(tab$log2FC) > logfc),]
  return(list(groupdata = group, resdf=tab, deg=deg))
}


#' Transformat DEG anlysis result obj to GSEA APP files
#' @description transformat DEG_edgeR obj to the files that GSEA APP need\cr
#' write the file that GSEA app need to local disk without anything retrun to object
#'
#' @param diffan.obj object DEG analysis result obj(such as dds)
#' @param outdir character the result output dir
#' @param startcol number the first data keep col, default is 8, please do not modify unless you know whate you are doing!
#' @param txt.del logical if delete the original txt file, default is TRUE, please do not modify unless you know whate you are doing!
#'
#' @return no
#' @export
#' @importFrom readr write_lines
#'
#' @examples
#' # suppressWarnings({
#' # DEG_resTogesa(diffan.obj = diffan.d,
#' #               outdir = z_result4)
#' # })
DEGres_ToGSEA <- function(diffan.obj, outdir, startcol = 8, txt.del=T) {
  outdir <- dirclean(outdir)
  mkdir.p(outdir)
  ### eset
  resdf <- diffan.obj$resdf
  df.gsea <- data.frame(NAME=resdf$Gene,
                        Description = 'Farmer',
                        resdf[, startcol:ncol(resdf)],
                        check.names = F)
  g.h <- c("#1.2" , paste(nrow(df.gsea), ncol(df.gsea)-2, sep = "\t"))
  readr::write_lines(g.h, file = paste0(outdir, "/gesa_eset", ".txt"),
              append = T)
  write.table(df.gsea, file = paste0(outdir, "/gesa_eset", ".txt"),
              sep = "\t", row.names = F, append = T)
  ### group
  group <- diffan.obj$groupdata[, 1, drop=F]
  g.g <- c(
    paste(nrow(group), length(unique(group[,1])), 1, sep = "\t"),
    paste0("#", "\t", paste(unique(group[,1]), collapse = "\t") ),
    paste(group[,1], collapse = "\t") )
  readr::write_lines(g.g, file = paste0(outdir, "/gesa_group", ".txt"))

  file.copy(paste0(outdir, "/gesa_eset", ".txt"),
            paste0(outdir, "/gesa_eset", ".gct"))
  file.copy(paste0(outdir, "/gesa_group", ".txt"),
            paste0(outdir, "/gesa_group", ".cls"))
  if (txt.del) {
    file.remove(paste0(outdir, "/gesa_eset", ".txt"))
    file.remove(paste0(outdir, "/gesa_group", ".txt"))
  }
}


#' Transformat DEG_edgeR obj TO res obj
#' @description transformat DEG_edgeR obj TO res obj for follow analysis
#'
#' @param diffan.obj DEG_edgeR obj
#' @param p pavlue
#' @param q fdr is q is set NULL, mean that fdr is not included for follow analysis
#' @param f logFC
#' @param mark outfile filename, if mark or outdir is not set then the result will not be written to local disk.
#' @param outdir outfile dirif mark or outdir is not set then the result will not be written to local disk.
#'
#' @return list(DIFF.ALL = res, DEG.ALL = genelist.deg,
#'              UP = genelist.up, DOWN = genelist.down,
#'              gsealist = genelist)
#'         DIFF.ALL: all deg table，genelist(deg,up,down)：genelist for enrichment
#'         gsealist：data for gsea analysis
#' @export
#' @importFrom writexl write_xlsx
#' @import dplyr
#'
#' @examples
#' #res <- DEG_edgeR.filter(diffobj = diffan.d, p = pvalue, q=fdr, f=log2(fc),
#' #                        mark = out_mark2, outdir = z_result2)
#' #
DEGres_ToRICH <- function(diffan.obj, p, q, f, mark, outdir) {
  # q如果等于NULL， 表示q不参与差异基因的筛选
  # vacalno 数据 【res】如果后两项参数缺失任意一个则只返回数据，不保存数据到硬盘
  res <- diffan.obj$resdf
  if(!is.null(q)) {
    deg <- res %>% filter(PValue < p & FDR < q & abs(log2FC) > f)
  } else {
    deg <- res %>% filter(PValue < p & abs(log2FC) > f)
  }
  res$Type <- "NONE"
  res[res$Gene %in% deg$Gene & res$log2FC >  f, "Type"] <- "UP"
  res[res$Gene %in% deg$Gene & res$log2FC < -f, "Type"] <- "DOWN"
  # gsea 数据 genelist
  res <- res %>% dplyr::arrange(desc(log2FC))
  genelist <- res$log2FC
  names(genelist) <- res$Gene
  # go kegg 数据 up down
  genelist.up <- res[res$Type == "UP", "Gene"] %>% as_tibble_col(column_name = "UP")
  genelist.down <- res[res$Type == "DOWN", "Gene"] %>% as_tibble_col(column_name = "DOWN")
  genelist.deg <- deg$Gene %>% as_tibble_col(column_name = "DEG.ALL")
  if( missing(mark) | missing(outdir) ) {
    return(list(DIFF.ALL = res, DEG.ALL = genelist.deg,
                UP = genelist.up, DOWN = genelist.down,
                gsealist = genelist))
  } else {
    outdir <- dirclean(outdir)
    mkdir.p(outdir)
    ## 保存结果
    tryCatch({
      write_xlsx(list(DIFF.ALL = res, DEG.ALL = genelist.deg,
                      UP = genelist.up, DOWN = genelist.down),
                 path = paste0(trimws(outdir, whitespace = "/"), "/DIFF.an_", mark, ".xlsx"))
    }, error = function(e) {
      print("保存结果至本地失败(不影响程序后继运行)")
    }, finally = {
      return(list(DIFF.ALL = res, DEG.ALL = genelist.deg,
                  UP = genelist.up, DOWN = genelist.down,
                  gsealist = genelist, diffan.obj = diffan.obj))
    })
  }
}


#' Title
#'
#' @param eset data.frame, the exprs data, the data must not be log2 transformed
#' @param group data.frame, the pheno information, the first coloumn is the type information and should be factor.
#' @param method character, the method of test, default is wilcox.test, alternative is t.test
#' @param pval number, the pvalue cutoff
#' @param fdr number, the fdr cutoff
#' @param logfc number, the log2(fc) value cutoff
#'
#' @return list(resdf = res, deg = deg)
#' @export
#' @importFrom qvalue qvalue
#' @importFrom dplyr filter arrange select `%>%`
#' @importFrom stats wilcox.test t.test
#' 
#' @examples
#' \dontrun{
#' eset <- data.table::fread("./protein_count_pn.csv", data.table = F)
#' eset <- eset %>% tibble::column_to_rownames(var = "Protein")
#' 
#' group <- read.csv("./group.csv", row.names = 1)
#' group$Type <- as.factor(group$Type)
#' 
#' diff <- DEG_tw.test(eset = eset, group = group)
#' }
#' @author Jiang
DEG_tw.test <- function(eset, group, method = "wilcox.test",
                        pval = 0.05, fdr = 0.1, logfc = log2(2)) {
  if (names(group)[1] != "Type") {
    warning("Just support the first coloumn and the first coloumn of group should be Type, the coloumn name will be changed to Type.")
    names(group)[1] <- "Type"
  }
  stopifnot(all(colnames(eset) == rownames(group)), is.factor(group$Type))
  stopifnot(method == "wilcox.test" | method == "t.test")

  if (method == "wilcox.test") {
    res <- apply(eset, 1, function(x) 
      wilcox.test(as.numeric(x) ~ as.character(group$Type))$p.value)
  } else {
    res <- apply(eset, 1, function(x) 
      t.test(as.numeric(x) ~ as.character(group$Type))$p.value)
  }
  #
  res <- cbind(PValue = res, eset)
  res$Gene <- rownames(res)
  duibi <- levels(group$Type)
  res[, duibi[1]] <- rowMeans(eset[, rownames(group)[group$Type == duibi[1]]])
  res[, duibi[2]] <- rowMeans(eset[, rownames(group)[group$Type == duibi[2]]])
  res[, duibi[1]] <- ifelse(res[, duibi[1]] == 0, 1, res[, duibi[1]])
  res[, duibi[2]] <- ifelse(res[, duibi[2]] == 0, 0.1, res[, duibi[2]])
  
  res$log2FC <- log2(res[, duibi[2]] / res[, duibi[1]])
  res$QValue <- qvalue::qvalue(res$PValue)$qvalues
  res <- res %>% dplyr::select(Gene, log2FC, PValue, QValue, 
                               duibi[1], duibi[2], everything()) %>% 
    dplyr::arrange(desc(log2FC), PValue)
  deg <- res %>% dplyr::filter(PValue < pval & QValue < fdr & abs(log2FC) > logfc)
  return(list(resdf = res, deg = deg))
}


