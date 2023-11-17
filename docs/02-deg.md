# RNAseq差异分析 {#deg}

## 差异分析 {#deg-mian}
### 加载包

```r
rm(list = ls());gc()
library(LZ)
library(tibble)
library(data.table)
library(DESeq2)
library(parallel)
library(BiocParallel)
library(ggplot2)
cat(" 您电脑线程为:", detectCores())
# 如果是12代以后的interCPU,建议最高不超过6或8。服务器可加大设置，但不能大于线程总数。
# 此处如果电脑性能一般，建议直接使用n=4或者2。
n = 4
# register(MulticoreParam(n)) 苹果和linux电脑使用这句替代下句
register(SnowParam(n))
```


### 数据准备
RNAseq下游分析必须准备两个文件：表达矩阵表格文件，分组表格文件<br>
将gene_count.csv，group.csv放在工作目录下<br>

* gene_count.csv 矩阵数据格式(数值型，整型)

| ID  | row1 | row2 | row3  | row4|
| :---  | :--- | :--- | :---  | :---|
| gene1 |  34  |  23  |  56   |  23 |
| gene2 |  35  |  23  |  12   |  23 |
| gene3 |  12  |  78  |  78   |  78 |

* group.csv 分组数据格式：需要组的行名 包含于 表达谱的列名 rownames(group) %in% colname(eset)

| Sample   |   Type    | BATCH|
| :----    | :----     | :----|
| rowname1 |   tumor   |  1   |
| rowname2 |   tumor   |  1   |
| rowname3 |   normal  |  2   |
| rowname4 |   normal  |  2   |

* 注意：<br>
表达文件中的基因名是SYMBOL还是ESembleID。如为EsembleID,要注意是有小数点的ID还是没有小数点的。
有小数点的形式为这样：ESEM00000123.34，没有点的是ESEM00000123。<br>
还有记住基因名这列的列名，建议统一设置为ID。<br>

### 文件夹准备
本文中有时也将文件夹称为目录，这两者等价。建议每个项目新建一个文件夹，例如本项目新建了一个名为
LZexample的文件夹，然后再在这个文件夹下建了一个data文件夹，以后data目录专门用来存原始文件，例
如RNAseq分析所需的eset.csv, group.csv或者更加原始的文件。<br>

**目录结构建议**：本项目的目录初始结构，建议每个项目按着这个形式来。项目文件夹LZexample这个文
件夹名建议改成有意义的名称，一眼便能看出这个项目是什么数据或者什么目的，而data文件夹名不建议更改。
图片


### 差异分析预设置

```r
# 设置工作目录，即之后所有的操作如果不指定文件夹，都将会在这个文件夹下进行
setwd("C:/data/LZexample")  # 按需更改成你的项目文件夹
getwd() # 检查是否更改工作目录成功了？
# Windows系统下默认的文件夹路径是 "C:/data/LZexample" 这中斜杆分隔文件夹，
# 如果是直接从win复制而来的，请将斜杠\改成反斜杠/, 就如下面设置的这样。（改成\\也行）

# 设置此次分析的标记
mark <- "T_C"  # 此次分析的标记(记录谁比谁或和筛选阈值，建议设置的有意义)
# 设置结果输出的文件夹，按需设定，可保持默认，如果修改只修改此处即可，下面无需修改。
# 第一次分析可以不用改，但如是第二分析，必须改这个，否则会覆盖第一次的结果。
outdir <- "result"  
outdirsub <- paste0(outdir, "/",mark)
outdirsub.gsea <- paste0(outdirsub, "/gsea")
outdirsub.rich <- paste0(outdirsub, "/rich")
# 差异分析阈值设定
ffdr <- 0.1
fpval <- 0.05
flogfc <- 1
```


### 差异分析

```r
# 读取并整理数据（如果都是按照上面的要求来的，不需要改这里的参数）
glist <- DEG_prepareData(eset_file = "data/eset.csv", # 表达数据的相对路径
                         group_file = "data/group.csv", # 分组文件的相对路径
                         id_dot = F,  # ESEM是否有点，有点设为T
                         col.by = "ID",  # 基因名列的列名
                         annot_trans = F, # 是否要注释，如果是EsembleID就需要设置为T
                         f_mark = mark)

# 差异分析 deseq2三部曲
dds <- DEG_DESeq2.dds(exprset.group=glist, batch = F)
DEG_DESeq2.pca(dds, outdir = outdirsub) # 此处有warning信息，不用管。
dds_list <- DEG_DESeq2.ana(dds)

# 差异后分析
if ( dir.exists(outdirsub.gsea) ) {
  # 如果outdirsub.gsea文件夹存在，清空该文件夹下所有文件
  file.remove(list.files(outdirsub.gsea, full.names = T))
}
# 构建GSEA官网软件分析所需格式文件
DEGres_ToGSEA(diffan.obj = dds_list, outdir = outdirsub.gsea) # 此处有warning，不用管。
# all_father中记录了
#             差异分析的总表，默认阈值的差异基因表，
#             上调基因列表，下调基因列表
#             以及R-GSEA分析所需要的所有mRNA的表达排序列表。
#   是我们后续各种分析的万恶之源（因此命名all_father）
#   上述数据同时被保存在本地硬盘的到一个多sheet的xlsx表中【outdirsub.rich目录中】
all_father <- DEGres_ToRICH(diffan.obj = dds_list, p=fpval, q=ffdr, f=flogfc,
                            mark=mark, outdir = outdirsub.rich)
save.image(file = paste0(outdirsub, "/1.diff.img.RDATA")) # 保存中间数据
```

## 火山图 {#deg-valcano}
需要修改的是以下几个值：<br/>
df_valcano：文件读取时的文件路径<br/>
ffdr: FDR阈值<br/>
fpval: PValue阈值<br/>
flogfc: logFC阈值<br/>
filterc: 火山图展示模式(p,fdr均考虑模式, 仅考虑p值模式，仅考虑fdr值模式，默认为第一种"fppadj")<br/>
label_gene: 展示基因列表<br/>
不清楚建议默认：fdr=0.1, pval=0.05, p,fdr均考虑模式。仅修改`label_gene: 展示基因列表`即可。

```r
# 本流程中不需要运行，后续想再次分析时可从此步开始
# load("./result/T_C/1.diff.img.RDATA") 
# library(LZ)
library(ggpubr);library(ggrepel);library(ggsci);library(scales)
library(tidyverse);library(dplyr);library(pheatmap);library(RColorBrewer)
# 导入火山图需要的数据，即差异分析后的未筛选表格(我们也称这个对象为resdf,
#  resdf文件涵盖差异分析的所有结果信息，可以做后续所有基于差异分析或者基因列表
#  的所有分析，如果后续分析时使用其它数据，请按这个resdf的格式改数据，主要
#  就是把数据的列名改成和resdf的列名相同，即可用此包的函数分析画图)
#  即df_valcano <- readxl::read_xlsx("xxx.xlsx", sheet = 1)
df_valcano <- all_father$DIFF.ALL
names(df_valcano) # 对应的列名必须为Gene, log2FC, PValue, FDR
# 差异分析阈值设定
ffdr <- 0.1
fpval <- 0.05
flogfc <- 1
# 模式设定
filterc <- "fppadj" # pvalue, padj均考虑模式 ("fpadj":仅考虑fdr值模式, "other": 仅考虑p值模式)
# 设定需要标记的marker gene
label_gene <- c('TFRC', 'ACSL1', 'LPCAT3', 'PCBP1', 'FTH1', 'SLC11A2',
                'SLC39A8', 'SAT1', 'FTL', 'GSS')
# 查看想展示的基因在不在差异分析总表中
# label_gene %in% df_valcano$Gene %>% all
# pic_data %>% filter(Row.names %in% label_gene)
# 火山图数据预处理
pic_data <- DEGp_prepareVolcano(df_valcano = df_valcano, filterc = filterc)
# 火山图 无标记
DEGp_Volcano(result = pic_data, logFC = flogfc,
             adj_P = ffdr, label_geneset = NULL)
ggsave(paste0(outdirsub,"/valcano.pdf"), width = 7, height = 7) # 保存
# 火山图 有标记
DEGp_Volcano(result = pic_data, logFC = flogfc, # log2(2)
             adj_P = ffdr, label_geneset = label_gene) %>%
  ggplotGrob() %>% cowplot::plot_grid()
ggsave(paste0(outdirsub, "/valcano.mark.gene.pdf"), width = 7, height = 7) # 保存
```

## 热图
On the way ...
