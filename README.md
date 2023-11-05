## LZ
本包致力于简化生信分析流程和批量分析。<br>
目前主要为RNAseq分析流程，后期会加入多组学联合分析流程。<br>

### 安装
BiocManager::install(c("DOSE", "clusterProfiler"))
devtools::install_github("ArronLZ/LZ", build_manual = T, build_vignettes = T, 
                         upgrade ="never", force = T)
