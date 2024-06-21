install_dependencies <- function() {
  cran_pack <- c('devtools', 'prettydoc', 'Hmisc', 'markdown', 
                 'Hmisc', 'tidyverse')
  for (p in cran_pack) { 
    if (!requireNamespace(p, quietly = T)) install.packages(p) 
  }
  
  bioc_pack <- c("DOSE", "GSVA", "biomaRt", "sva", "quantiseqr",
                 "singscore", "preprocessCore", 'DESeq2', 'edgeR', 
                 'limma', "topGO", 'Rgraphviz', 'org.Hs.eg.db', 
                 "clusterProfiler")
  for (p in bioc_pack) { 
    cat(p, '=========\n') 
    if (!requireNamespace(p, quietly = T)) 
      BiocManager::install(p, update = F, ask =F) 
  }
}

install_dependencies()