## ----eval=FALSE---------------------------------------------------------------
#  rm(list = ls());gc()
#  library(LZ)
#  library(tibble)
#  library(parallel)
#  library(data.table)
#  library(DESeq2)
#  cat(" 您电脑线程为:", detectCores())
#  # 数字建议设置为线程数的一半，高配电脑可以加大设置，但不能大于电脑线程总数
#  n = round(detectCores()/2)  # 此处如不懂不要修改。
#  # register(MulticoreParam(n)) 苹果和linux电脑使用这句替代下句
#  register(SnowParam(n))

## ----eval=FALSE---------------------------------------------------------------
#  outdir <- "result"  # 按需设定（可保持默认，如果修改只修改此处即可，下面无需修改）
#  outdirsub <- paste0(outdir, "/1.diff")
#  outdirsub.gsea <- paste0(outdirsub, "/gsea")
#  outdirsub.rich <- paste0(outdirsub, "/rich")

## ----eval=FALSE---------------------------------------------------------------
#  glist <- DEG_prepareData(eset_file = "gene_count.csv",
#                           group_file = "group.csv",
#                           annot_trans = T,
#                           f_mark = "SHUANG-con")
#  # 差异分析 deseq2三部曲
#  dds <- DEG_DESeq2.dds(exprset.group=glist, batch = F)
#  DEG_DESeq2.pca(dds, outdir = outdirsub) # 此处有warning信息，不用管。
#  dds_list <- DEG_DESeq2.ana(dds)
#  
#  # 差异后分析
#  if ( dir.exists(outdirsub.gsea) ) {
#    file.remove(list.files(outdirsub.gsea, full.names = T))
#  }
#  DEGres_ToGSEA(diffan.obj = dds_list, outdir = outdirsub.gsea) # 此处有warning信息，不用管。
#  all_father <- DEGres_ToRICH(diffan.obj = dds_list, p=0.05, q=0.1, f=1,
#                              mark="SHUANG-CONTROL", outdir = outdirsub.rich)
#  save.image(file = "./result/1.diff.img.RDATA")

## ----eval=FALSE---------------------------------------------------------------
#  # rm(list = ls());gc()
#  #library(ggpubr);library(ggrepel);library(ggsci);library(scales)
#  library(tidyverse);library(dplyr);library(pheatmap);library(RColorBrewer)
#  library(xlsx)
#  # # 导入火山图需要的差异分析后的基因全部表格
#  df_valcano <- xlsx::read.xlsx("./result/1.diff/rich/DIFF.an_SHUANG-CONTROL.xlsx", sheetIndex = 1)
#  # saveRDS(df_valcano, file = "./result/1.diff/rich/DIFF.an_SHUANG-CONTROL.rds")
#  # df_valcano <- readRDS("./result/1.diff/rich/DIFF.an_SHUANG-CONTROL.rds")
#  names(df_valcano) # 对应的列名必须为Gene, log2FC, PValue, FDR
#  # 阈值设定
#  ffdr <- 0.2
#  fpval <- 0.05
#  flogfc <- 1
#  # 模式设定
#  filterc <- "fppadj" # pvalue, padj均考虑模式 ("fpadj":仅考虑fdr值模式, "other": 仅考虑p值模式)
#  # 火山图数据预处理
#  pic_data <- DEG_prepareVolcano(df_valcano = df_valcano, filterc = filterc)
#  
#  # 设定需要标记的marker gene
#  label_gene <- c('TFRC', 'ACSL1', 'LPCAT3', 'PCBP1', 'FTH1', 'SLC11A2',
#                  'SLC39A8', 'SAT1', 'FTL', 'GSS')
#  # 查看想展示的基因在不在差异分析总表中
#  # label_gene %in% df_valcano$Gene %>% all
#  # pic_data %>% filter(Row.names %in% label_gene)
#  
#  # 火山图 无标记
#  DEGplot_Volcano(result = pic_data, logFC = flogfc,
#                 adj_P = ffdr, label_geneset = NULL)
#  # 保存图片
#  # ggsave("./valcano.pdf", width = 7, height = 7)
#  # 火山图 有标记
#  DEGplot_Volcano(result = pic_data, logFC = flogfc, # log2(2)
#                 adj_P = ffdr, label_geneset = label_gene) %>%
#    ggplotGrob() %>% cowplot::plot_grid()
#  # 保存图片
#  # ggsave("./result/1.diff/valcano.mark.gene.pdf", width = 7, height = 7)

