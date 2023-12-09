#' @title S3 DEG_prepareData: general function DEG_prepareData
#' @description S3 method initialize, initialize the S3 general function DEG_prepareData()
#' 
#' @param obj object
#'
#' @return list
#' @export
DEG_prepareData <- function(obj) { UseMethod("DEG_prepareData", obj) }

#' @title S3 DEG_prepareData.default: prepare the data that LZ diff analysis need
#' @description S3 method of default, prepare the data that LZ diff analysis need, \cr
#' #' return a list(eset=xx, group=xx, f_mark="xx"), S3 general function DEG_prepareData.default
#' @param eset_file character, the exps data file name, default: "gene_count.csv"
#' @param id_dot logical, if the id column is ensemble id with dot, ex.ESEMxxxx.3
#' @param col.by character, delete the duplicate by this column
#' @param col.del character, manually delete some column
#' @param auto.del.character logical, if auto delete the column those class is character
#' @param group_file the group data file name, default: "group.csv"
#' @param annot_trans logical, do the eset need be annot to symbol? if the datafame's first coloumn is not offical symbol, you can set to TRUE.
#' @param f_mark file mark(optional), default: diff
#' @param is.human logical, default is True, if it is True, the param orgDB,fromTyp and toType will be ignore
#' @param orgDb character, annotation db, ex. "org.Mm.eg.db"
#' @param fromType character, eset id column type
#' @param toType character, the annot result
#'
#' @return list
#' @export
#' @importFrom data.table fread
#' @importFrom clusterProfiler bitr
#'
#' @examples #
#' #' # exprset.group = list(eset, group, f_mark)
#' # ** 如果bacth=T，则group中第二列为batch列。
#' # ** 矩阵数据格式(数值型，整型)
#' # gene  | row1 | row2 | row3  | row4
#' # gene1 |  34  |  23  |  56   |  23
#' # gene2 |  35  |  23  |  12   |  23
#' # gene3 |  12  |  78  |  78   |  78
#' # ** 分组数据格式：需要组的行名 包含于 表达谱的列名 rownames(group) %in% colname(eset)
#' # Sample   |   Type    | BATCH
#' # rowname1 |   tumor   |  1
#' # rowname2 |   tumor   |  1
#' # rowname3 |   normal  |  2
#' # rowname4 |   normal  |  2
#' # if is.human = F, please set the orgDb and fromType correctly.
DEG_prepareData.default <- function(eset_file="gene_count.csv",
                            id_dot = F, col.by = "ID",
                            col.del=NULL, auto.del.character=T,
                            group_file="group.csv",
                            annot_trans=T, f_mark="diff",
                            is.human = T, orgDb = "org.Mm.eg.db",
                            fromType = "ENSEMBL", 
                            toType = c("SYMBOL", "UNIPROT")
                            ) {
  eset <- data.table::fread(eset_file, data.table = F)
  eset <- quchong(eset = eset, col.by = col.by, col.del = col.del, 
                  auto.del.character = auto.del.character)
  
  # if the id contain ., only keep the part before .  ---------
  if (id_dot) {
    eset$annot <- sapply(strsplit(rownames(eset), "\\."), 
                         function(x) x[[1]])
    eset <- quchong(eset = eset, col.by = "annot")
  }
  # if the id contain ., only keep the part before .///  ---------
  
  # id annot or not  ----------
  if (is.human) {
    # human
    tmpenv <- new.env()
    data(gencode.v22.annot, package = "LZ", envir = tmpenv)
    all_anot <- tmpenv$gencode.v22.annot
    if (annot_trans) {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$ensembl_gene_id, 
      ]
      eset <- cbind(SYMBOL = all_anot$hugo_mRNA[match(rownames(eset), 
                                                      all_anot$hugo_mRNA$ensembl_gene_id), 3], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    } else {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$symbol, 
      ]
    }
    rm(tmpenv)
  } else {
    # not human , default is mm
    if (annot_trans) {
      annot.df.org <- bitr(rownames(eset), fromType = fromType, 
                           toType = toType, OrgDb = orgDb)
      eset <- eset[rownames(eset) %in% annot.df.org[, 1], ]
      eset <- cbind(SYMBOL = annot.df.org[match(rownames(eset), 
                                                annot.df.org[, 1]), "SYMBOL"], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    }
  }
  # id annot or not ///  ----------
  
  group <- read.csv(group_file, row.names = 1)
  
  # check if the rowname of group is all from the colnames of eset -----
  if ( all(rownames(group) %in% names(eset)) ) {
    # eset's coloumns > group's rows
    if (length(names(eset)) > length(rownames(group))) {
      cat(" 请注意：!!!\n 提供的分组文件中样本数量少于表达矩阵中的样本数：\n
 这意味着将会只保留分组文件中的样本来进行后续分析，\n 不在分组文件中的样本将会被剔除!!!\n
 分组文件内容展示：\n")
      print(group)
    }
    # only keep the eset's coloumn by group's rownames sample
    eset <- eset[, rownames(group)]
    # check all eset's coloumn == group's rowname
    if (all(names(eset) == rownames(group))) {
      glist <- list(eset = eset, group = group, f_mark = f_mark)
      glist$group[,1] <- factor(glist$group[,1])
    } else {
      stop("提供的分组信息和表达信息不匹配，请检查后再运行\n")
    }
  } else {
    stop("提供的分组信息和表达信息不匹配(分组文件有样本名不在表达矩阵的样本名中)，\n
         请检查后再运行\n")
  }
  # check if the rowname of group is all from the colnames of eset /// -----
  return(glist)
}

#' @title S3 DEG_prepareData.DEres: prepare the data that LZ diff analysis need
#' @description S3 method of object DEres, prepare the data that LZ diff analysis need, \cr
#' return a list(eset=xx, group=xx, f_mark="xx"), S3 general function DEG_prepareData.DEres
#' 
#' @param obj object
#' @param f_mark character a analysis mark(optional), default: "DE"
#'
#' @return list
#' @export
DEG_prepareData.DEres <- function(obj, f_mark = "DE") {
  glist <- list(eset = obj$eset, group = obj$group, f_mark = f_mark)
  return(glist)
}


#' RNAseq differential analysis DESeq2 001<dds>
#' @description The function is designed to perform differential analysis by DESeq2<dds>.
#' @param exprset.group the list consisted of eset, group and mark.
#' @param batch batch number
#'
#' @return dds ( DESeq(dds, parallel = T) )
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @import BiocParallel
#'
#' @examples
#' # exprset.group = list(eset, group, f_mark)
#' # ** 如果bacth=T，则group中第二列为batch列。
#' # ** 矩阵数据格式(数值型，整型)
#' #       | row1 | row2 | row3  | row4
#' # gene1 |  34  |  23  |  56   |  23
#' # gene2 |  35  |  23  |  12   |  23
#' # gene3 |  12  |  78  |  78   |  78
#' # ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eset)
#' #            Type   BATCH
#' # rowname1 | tumor  1
#' # rowname2 | tumor  1
#' # rowname3 | normal 2
#' # rowname4 | normal 2
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
  dds <- DESeqDataSetFromMatrix(countData = exprset,
                                colData = pheno,
                                design = design_type)
  keep <- rowSums(counts(dds) >= 10) >= ncol(exprset)/3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数
  dds <- dds[keep, ]
  # 1.3
  dds <- DESeq(dds, parallel = T)
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
#' @importFrom DESeq2 vst plotPCA
#' @importFrom factoextra hcut fviz_dend
#' @importFrom stats dist
#'
#' @examples
#' #
DEG_DESeq2.pca <- function(dds, outdir) {
  # if (!dir.exists(outdir)) {
  #   dir.create(outdir, recursive = T, showWarnings = T)
  # }
  outdir <- dirclean(outdir)
  mkdir.p(outdir)
  ## 数据转换
  # library(factoextra)
  vsd <- vst(dds, blind=FALSE)
  plotpca <- plotPCA(vsd, intgroup="Type")
  sampleDists <- dist(t(assay(vsd)))
  h.data <- hcut(sampleDists, k = 2, stand = FALSE, hc_method = "average")
  ### plot
  pdf(paste0(outdir, '/pre_PCA&Dispersion.plot', '.pdf'), width = 8, height = 7)
  # 1. Dispersion plot
  plotDispEsts(dds, main="Dispersion Plot")
  pf <- fviz_dend(h.data, rect_fill = T, cex = 0.8, color_labels_by_k = T, horiz = T)
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
  plotDispEsts(dds, main="Dispersion Plot")
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
#' @import dplyr
#'
#' @examples
#' #
DEG_DESeq2.ana <- function(dds, pval=0.05, fdr=0.1, logfc=1) {
  # resultsNames(dds)                # 查看结果的名称。
  # dds$type                         # 默认后比前
  type1 <- levels(dds$Type)[1]
  type2 <- levels(dds$Type)[2]
  # 如不指定contrast，默认第二因子比第一因子。建议将对照的因子水平放在第一个。
  # 也可以通过results(dds, contrast=c("type", type2, type1)): type2 比 type 1
  res <- results(dds)
  # summary(res)      #看一下结果的概要信息，p值默认小于0.1。

  ### 差异分析总表
  res_df <- data.frame(res)
  eset_norma <- data.frame(counts(dds, normalized=TRUE), check.names = F)
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
#' @importFrom stats median model.matrix
#' @import edgeR
#' @import dplyr
#'
#' @examples # Run: diffan <- DEG_edgeR(exprset.group) # pval=0.05, fdr=0.1, logfc=1
#' # exprset.group  由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
#' # *1 exprset.group$eset： 矩阵数据格式(数值型，整型)
#' #       | row1 | row2 | row3  | row4
#' # gene1 |  34  |  23  |  56   |  23
#' # gene2 |  35  |  23  |  12   |  23
#' # gene3 |  12  |  78  |  78   |  78
#' # *2 exprset.group$group： 分组数据格式：建议只有一列数据
#' #    # 需要组的行名=表达谱的列名 rownames(group) == colname(eset)
#' #            type
#' # rowname1 | tumor
#' # rowname2 | tumor
#' # rowname3 | normal
#' # rowname4 | normal
#' # *3 exprset.group$f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
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
  y <- DGEList(counts = exprset, group = pheno[,1])
  keep <- filterByExpr(y)
  #keep <- filterByExpr(y, min.count = 10, min.total.count = 15,
  #                     large.n = 10, min.prop = 0.7)
  y <- y[keep, , keep.lib.sizes=FALSE]
  # 进行TMM标准化
  y.norm <- calcNormFactors(y, method = "TMM")
  # 获取cpm TMM值
  cpm.tmm <- cpm(y.norm, log = F)
  cpm.tmm <- log2(cpm.tmm + 1)

  # 1.2 计算离散因子
  y <- calcNormFactors(y)
  # design.form <- as.formula(paste("~", "group"))
  # # y$samples$group <- relevel(y$samples$group, ref="C")
  design <- model.matrix(~group, data =  y$samples)
  y <- estimateDisp(y, design)
  #To perform quasi-likelihood F-tests:
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef)

  # 获取结果
  # 获取排名靠前的基因，这里设置n=80000是为了输出所以基因
  et <- topTags(qlf, n = 80000)
  et <- as.data.frame(et) # 转换为数据框类型
  et.norm <- merge(et, cpm.tmm, by = 0)
  et.norm <- column_to_rownames(et.norm, var = "Row.names")
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
#'
#' @examples
#' #
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
  y <- calcNormFactors(y)
  # voom
  v <- voom(y, design, plot=F)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  #topTable(fit, coef=ncol(design))
  tab <- topTable(fit, sort.by = "P", n = Inf)
  tab <- cbind(rownames(tab), tab)  # 将行名粘贴为数据框的第一列
  colnames(tab) <- c("Gene", "log2FC", "AveExp", "T", "PValue", "FDR", "B")
  tab <- tab %>% dplyr::select(Gene, log2FC, PValue, FDR, everything())
  tab <- tab %>% dplyr::arrange(desc(log2FC), PValue)

  deg <- tab[which(tab$PValue < pval & tab$FDR < fdr & abs(tab$log2FC) > logfc),]
  return(list(groupdata = pheno, resdf=tab, deg=deg))
}


#' General compare the differnt ID(Gene) between the group B vs A
#' @description General compare the differnt ID(Gene) between the group B vs A
#'
#' @param eset data.frame, the exprs data
#' @param group data.frame, the pheno information, the first coloumn is the type information and should be factor.
#'
#' @return # all deg dataframe(resdf)
#' @export
#' @importFrom stats model.matrix
#' @importFrom limma lmFit eBayes topTable
#' @import dplyr
#'
#' @examples #
#' # pn_eset <- data.table::fread("./protein_count_pn.csv", data.table = F)
#' # pn_eset <- pn_eset %>% tibble::column_to_rownames(var = "Protein")
#' #
#' # group <- read.csv("./group.csv", row.names = 1)
#' # group$Type <- as.factor(group$Type)
#' #
#' # df <- log2(pn_eset + 1)
#' # diff <- limma.general(eset = df, group = group)
limma.general <- function(eset, group) {
  names(group)[1] <- "Type"
  group_list <- group$Type
  design <- model.matrix(~group_list) # ~法,截距法
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(eset)
  # ~三部曲
  fit <- lmFit(eset, design)
  fit <- eBayes(fit, trend=TRUE)
  tab <- topTable(fit, coef=2, n=Inf)
  tab <- merge(tab, eset, by = 0)
  tab <- tab %>%
    dplyr::rename(ID=Row.names, log2FC=logFC, PValue=P.Value, FDR=adj.P.Val) %>%
    dplyr::select(ID, log2FC, PValue, FDR, names(eset), everything()) %>%
    dplyr::arrange(desc(log2FC), PValue)
  return(tab)
}


#' Transformat DEG_edgeR obj to GSEA APP files
#' @description transformat DEG_edgeR obj to the files that GSEA APP need
#'
#' @param diffan.obj DEG_edgeR obj
#' @param outdir the result output dir
#'
#' @return write to local disk without anything retrun to object
#' @export
#' @importFrom readr write_lines
#'
#' @examples
#' # suppressWarnings({
#' # DEG_resTogesa(diffan.obj = diffan.d,
#' #               outdir = z_result4)
#' # })
DEGres_ToGSEA <- function(diffan.obj, outdir) {
  # diffan.obj # DEG_edgerR结果，outdir：结果输出文件夹，不能带/。
  # if (!dir.exists(outdir)) {
  #   dir.create(outdir, recursive = T, showWarnings = T)
  # }
  outdir <- dirclean(outdir)
  mkdir.p(outdir)
  ### eset
  resdf <- diffan.obj$resdf
  df.gsea <- data.frame(NAME=resdf$Gene,
                        Description = 'Farmer',
                        resdf[, 7:ncol(resdf)],
                        check.names = F)
  g.h <- c("#1.2" , paste(nrow(df.gsea), ncol(df.gsea)-2, sep = "\t"))
  write_lines(g.h, file = paste0(outdir, "/gesa_eset", ".txt"),
              append = T)
  write.table(df.gsea, file = paste0(outdir, "/gesa_eset", ".txt"),
              sep = "\t", row.names = F, append = T)
  ### group
  group <- diffan.obj$groupdata[, 1, drop=F]
  g.g <- c(
    paste(nrow(group), length(unique(group$group)), 1, sep = "\t"),
    paste0("#", "\t", paste(unique(group$group), collapse = "\t") ),
    paste(group$group, collapse = "\t") )
  write_lines(g.g, file = paste0(outdir, "/gesa_group", ".txt"))
  file.copy(paste0(outdir, "/gesa_eset", ".txt"),
            paste0(outdir, "/gesa_eset", ".gct"))
  file.copy(paste0(outdir, "/gesa_group", ".txt"),
            paste0(outdir, "/gesa_group", ".cls"))
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
                  gsealist = genelist))
    })
  }
}





