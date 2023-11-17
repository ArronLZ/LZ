# 安装 {#install}

**LZ** R包可以从Github上安装。<br>
先安装R及Rstudio, 前者是核心，后者是编辑器(写代码的地方)。<br>
1. 安装[R最新版](https://www.r-project.org/)<br>
2. 安装[Rstudio最新版](https://posit.co/download/rstudio-desktop/)<br>
3. Win电脑可以考虑安装R版本对应的[Rtools](https://cran.r-project.org/bin/windows/Rtools/) (可选项)

```r
# 安装biocondutor
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  # R 4.2.0
  BiocManager::install(version = "3.18") # 4.3 == 3.18
}

# 设置镜像(可选)
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 安装LZ依赖包
# 安装cran包
cran_pack <- c('devtools', 'prettydoc', 'Hmisc', 'markdown', 'tidyverse')
for (p in cran_pack) { if (!requireNamespace(p, quietly = T)) install.packages(p) }
# 安装bioconductor包
bioc_pack <- c("HPO.db", "DOSE", "clusterProfiler", 'DESeq2', 'edgeR', 'limma')
for (p in bioc_pack) { 
  cat(p, '=========\n') 
  if (!requireNamespace(p, quietly = T)) BiocManager::install(p, update = F, ask =F) 
}
# 安装LZ包
devtools::install_github("ArronLZ/LZ", upgrade ="never", force = T, 
                         build_vignettes = T)
# 查看文档（点击LZ Documents查看网页版）
vignette('LZ')
```
