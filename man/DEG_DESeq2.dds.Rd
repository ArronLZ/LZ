% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.R
\name{DEG_DESeq2.dds}
\alias{DEG_DESeq2.dds}
\title{RNAseq differential analysis DESeq2 001<dds>}
\usage{
DEG_DESeq2.dds(exprset.group, batch = F)
}
\arguments{
\item{exprset.group}{the list consisted of eset, group and mark.}

\item{batch}{batch number}
}
\value{
dds ( DESeq(dds, parallel = T) )
}
\description{
The function is designed to perform differential analysis by DESeq2<dds>.
}
\examples{
# Run: diffan <- DEG_DESeq2.dds(exprset.group) # batch = F
# exprset.group 由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
# ** 如果bacth=T，则group中第二列为batch列。
# 1. eset：
# ** 表达谱数据格式(数值型，整型)
#       | sample1 | sample2 | sample3  | sample4
# gene1 |  34     |  23     |  56      |  23
# gene2 |  35     |  23     |  12      |  23
# gene3 |  12     |  78     |  78      |  78
# 2. group：
# ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eseti)
# BATCH列可以没有(如果设置了bacth=T，BATCH列必须有，且必须位于第二列)
#            Type     BATCH
# sample1  | tumor   | 1
# sample2  | tumor   | 1
# sample3  | control | 2
# sample4  | control | 2
# 3. f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
# 默认的比法是后比前（字符排序在后面的比前面的，如示例中的就是 t vs c,因为c排在前面）

# LZ包将DESeq2差异分析整合拆解成了三个函数，因此需要如下代码才能完整进行差异分析。
# 差异分析 deseq2三部曲
# dds <- DEG_DESeq2.dds(exprset.group=glist, batch = F)  # 差异分析1：准备工作
# DEG_DESeq2.pca(dds, outdir = outdirsub)   # 画PCA图。 此处有warning信息，不用管。
# dds_list <- DEG_DESeq2.ana(dds)           # 差异分析1：得到结果并进行初步筛选
}
\author{
Jiang
}
