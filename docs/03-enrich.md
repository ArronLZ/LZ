# 富集分析 {#enrich}

## GO & KEGG 富集分析 {#enrich-auto}
### 一键脚本(批量处理)
这是一个一键脚本，请新建一个单独的文件写这段脚本，然后按这个脚本的顶部注释修改`resdf outd fc.list处即可`，运行即可批量出不同FC的富集分析结果。

```r
# 此脚本为GO、KEGG分析（需要一个输入文件即可，为差异分析流程后的resdf文件）
# 即为第一步(或1脚本)的结果的一个结果文件（DIFF_an_***.xlsx）
# 即为resdf文件，此文件是差异分析后的总表
# 注意如果采用了其他的分析方法得到差异分析后表，运行这个脚本时可能需要更改列名
# 即我们的resdf对象的列名为Gene, log2FC,PValue,FDR，需要与这些个列名保持一致。
# 此脚本中的需要修改的位于 /// *** /// 行中，另外还有一个LZ::setproxy()行，
#   如果没有代理工具，或者代理工具不支持http代理，或者端口不通，请不要运行。
rm(list = ls());gc() # 清空所有对象，慎用，必要时用
suppressMessages({ suppressWarnings({
  library(LZ)
  library(tidyverse);library(data.table)
  library(clusterProfiler);library(enrichplot)
  library(topGO);library(Rgraphviz)
  library(RColorBrewer);library(ggsci);library(pheatmap)
  library(xlsx);library(readxl)
}) })
# 若无代理工具，切勿运行 
# LZ::setproxy() # 高危！！！新手不要运行此行，会使当前窗口断网！！！
# Sys.getenv('http_proxy') Sys.setenv('http_proxy'='') Sys.setenv('https_proxy'='')



# /// 1,2,3均需按实际改写
# 读取数据  resdf存放目录
resdf <- readxl::read_xlsx("result/rnaseqOR-NC/rich/DIFF.an_OR-NC.xlsx",
                  sheet = 1) %>% as.data.frame()
# 此处可能需要插入修改列名的代码，需要为标准的resdf格式
#   标准resdf格式，用列名Gene, log2FC, PValue, FDR来表示gene, log2fc, p, q/fdr

# 输出目录
outd = "result/rnaseqOR-NC/rich" 
# logFC阈值, 多个阈值的话，写成fc.list <- list('1.2' = log2(1.2), '2' = log2(2))
# 注意！！！！！！：括号里log2(2)的2，和引号里'2'的2都要需同步要改。！！！
# 否则可能会覆盖结果
fc.list <- list('2'=log2(2))
# 设置物种为人类（如是人类则不需要更改）
GO_database <- 'org.Hs.eg.db' # keytypes(org.Hs.eg.db)
KEGG_database <- 'hsa' 
# ///


# 预处理数据符合GOKEGG分析的要求
# # 不同fc条件下的GOgenelist list(ALL, UP, DOWN)
gogenelist <- lapply(fc.list, function(x) DEG_prepareGOglist(resdf, logfc = x))
# gogenelist %>% length()
#
# 对logFC迭代，每个FC新建一个目录，用来存upgene, downgene, allgene的GO结果
go <- DEG_runGO(outdir = outd, genelist = gogenelist)
# 对logFC迭代，每个FC新建一个目录，用来存upgene, downgene, allgene的KEGG结果
kegg <- DEG_runKEGG(outdir = outd, genelist = gogenelist)
```

### 简易GO,KEGG一次分析  {#enrich-simple}
如果已经得到了差异基因列表，且无需批量分析，可以进行这个简易分析。<br/>
数据格式: `head(genelist.lh)` <br/>
[1] "AARS1"  "AATF"   "ABCB7"  "ABCE1"  "ABHD11" "ABHD12"

```r
# 简易GO,KEGG一次分析(即：已经得到了差异基因列表)
# LZ::setproxy() # 代理设置，新手别碰，会断网
# 差异基因列表
genelist.lh <- pic.list$sig.data$Gene 
# 转换ID
gene_df <- bitr(genelist.lh, fromType = "SYMBOL", toType = c("ENTREZID", "UNIPROT"), 
                OrgDb = 'org.Hs.eg.db')
# GO分析
go.lh <- DEG_GO(gene_df, orgdb = "org.Hs.eg.db", sigNodes = 20, 
                resultdir="./result/proteinOR-NC", filemark = "p1.5_g_2")
go.lhdf <- sapply(go.lh, function(x) x@result, simplify = T)
write_xlsx(go.lhdf, path = "./result/proteinOR-NC/lh_go.all.xlsx")
# KEGG分析
kegg.lh <- DEG_KEGG(gene_df)
write_xlsx(kegg.lh$pSigDF, path = "./result/proteinOR-NC/lh_kegg.all.xlsx")
```

### GO、KEGG分析结果可视化 {#enrich-visual} {#dotplot}

```r
# dotplot go
# 读取go分析保存的表格
dotData <- readxl::read_xlsx("./result/proteinOR-NC/lh_go.all.xlsx", sheet = 1)
# 筛选数据（按需配合其他筛选）
dotData <- DEGp_prepareDotplot(dotData, head = 30, delete = NULL)
pic.dot <- DEGp_Dotplot(dotData, title = 'TOP of GO', 
                        resultdir = "./result/proteinOR-NC", filemark = 'GO_top', 
                        pic.save = T)

# dotplot kegg
# 读取kegg分析保存的表格
dotDatak <- readxl::read_xlsx("./result/proteinOR-NC/lh_kegg.all.xlsx", sheet = 1)
# 筛选数据（按需配合其他筛选）
dotDataK <- DEGp_prepareDotplot(dotDatak, head = 30, delete = NULL)
pic.dotk <- DEGp_Dotplot(dotDataK, title = 'TOP of KEGGpathway', 
                         resultdir = "./result/proteinOR-NC", filemark = 'KEGG_top', 
                         pic.save = F)

# 组合图
gh <- ggplotGrob(pic.dot)
gd <- ggplotGrob(pic.dotk)
cowplot::plot_grid(gh, gd, rel_widths = c(1, 1.25))
ggsave(paste0(dir_out, "/GO_KEGG_top.pdf"), width = 16, height = 10)
```

## GSEA分析 {#enrich-gsea}
### R GSEA
On the way ...

### GSEA官方软件
On the way ...
