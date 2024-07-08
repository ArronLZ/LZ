## LZ
本包致力于简化生信分析流程和批量分析。<br>
目前主要为RNAseq分析，后期会加入其他分析流程。<br>

**安装**
```r
BiocManager::install(c("DOSE", "clusterProfiler"))
devtools::install_github("ArronLZ/LZ", upgrade ="never", build_vignettes = F, force = T)
```

**文档** <br>
[点我进入说明文档](https://arronlz.github.io/LZ/)
目前完有RNAseq和ATAC分析流程 <br>
