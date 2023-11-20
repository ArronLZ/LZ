# 安装 {#install}

## 安装R及Rsudio环境
**LZ** R包可以从Github上安装。<br>
先安装R及Rstudio, 前者是核心，后者是编辑器(写代码的地方)。<br>
1. 安装[R最新版](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/)，根据系统自行选择
版本，win用户可以直接点[R-4.3.2](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/R-4.3.2-win.exe)下载。<br>
2. 安装[Rstudio最新版](https://posit.co/download/rstudio-desktop/)，win用户可直接点[Rstudio Desktop](https://download1.rstudio.org/electron/windows/RStudio-2023.09.1-494.exe)下载<br>
3. Win电脑可以考虑安装R版本对应的[Rtools](https://cran.r-project.org/bin/windows/Rtools/) (可选项，新手可以不安装)<br>
4. **重要提示：请卸载或者至少退出一切杀毒软件(微软自带的不用退出)，否则安装包时可能会出现难以解决的奇怪bug。**

## 安装LZ包
1. 安装LZ包所需的依赖包

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
cran_pack <- c('devtools', 'prettydoc', 'Hmisc', 'markdown', 
               'Hmisc', 'tidyverse')
for (p in cran_pack) { if (!requireNamespace(p, quietly = T)) install.packages(p) }
# 安装bioconductor包
bioc_pack <- c("HPO.db", "DOSE", "clusterProfiler", 'DESeq2', 'edgeR', 
               'limma', "topGO", 'Rgraphviz', 'org.Hs.eg.db')
for (p in bioc_pack) { 
  cat(p, '=========\n') 
  if (!requireNamespace(p, quietly = T)) BiocManager::install(p, update = F, ask =F) 
}
```

2. 安装LZ包(此包会不定时更新，后续更新只需要重新运行这句即可，上面的包不需要重新安装)

```r
# 安装LZ包
devtools::install_github("ArronLZ/LZ", upgrade ="never", force = T, 
                         build_vignettes = T)
# 查看文档（点击LZ Documents查看网页版）
vignette('LZ')
```
