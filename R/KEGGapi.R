#' Get KEGG compound
#' @description 获取指定人类hsa通路中的代谢物（通过R KEGGREST实现）
#' @param id character(1), hsa_id such as "hsa00010", ...
#' 
#' @return dataframe
#' @export
#' @importFrom KEGGREST keggGet
#' @author Jiang
getKEGGcpd <- function(id) {
  # id = names(pathway)[1]
  # id = "hsa00010"
  hsa <- KEGGREST::keggGet(id)
  compound <- hsa[[1]]$COMPOUND
  path.class <- ifelse(is.null(hsa[[1]]$CLASS), NA, hsa[[1]]$CLASS)
  path.name <- ifelse(is.null(hsa[[1]]$NAME), NA, hsa[[1]]$NAME)
  if (!is.null(compound)) {
    cpd.name <- as.character(compound)
    cpd.id <- names(compound)
    metabolism.df <- data.frame(
      ID = id,
      NAME = path.name,
      CLASS = path.class,
      cpd.id = cpd.id,
      cpd.name = cpd.name)
  } else {
    cat(id, "have no compound\n")
    return(NULL)
  }
  return(metabolism.df)
}

#' Get KEGG compound Loop
#' @description 获取所有人类hsa通路中的代谢物(通过循环实现)
#' @param pathway character(n), name character the name of `pathway` is the hsa_id
#'
#' @return dataframe
#' @export
#' @importFrom stringr str_split
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr `%>%`
#' @author Jiang
getKEGGcpd.loop <- function(pathway) {
  hsa_id <- names(pathway)
  re <- list()
  for (i in seq_along(hsa_id)) {
    if((i %% 100) == 0) {
      cat("Stop 10 second...\n")
      Sys.sleep(10)
    } 
    cat(i, hsa_id[i], "\n")
    df <- getKEGGcpd(id = hsa_id[i])
    if (!is.null(df)) {
      re[[hsa_id[i]]] <- df
    }
    Sys.sleep(1)
  }
  
  re.df <- Reduce(rbind, re)
  re.df <- cbind(re.df, str_split(re.df$CLASS, ";", simplify = T))
  re.df$CLASS <- NULL
  # 
  re.df <- re.df %>% 
    dplyr::rename(category = "1", subcategory = "2") %>% 
    dplyr::select(category, subcategory, everything())
  return(re.df)
}


#' Get KEGG compound
#' @description 获取所有代谢通路及代谢物对应关系(通过原生kegg api实现)
#' @return list
#' @export
#' @importFrom readr read_tsv
#' @importFrom dplyr `%>%`
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @author Jiang
getKEGG.metabolism <- function() {
  # aqurie data from online 
  tryCatch({
    gmt.path.cpd <- read_tsv("https://rest.kegg.jp/link/cpd/pathway", col_names = F) %>% 
      data.frame(check.names = F)
    cpd.meta <- read_tsv('https://rest.kegg.jp/list/compound', col_names = F) %>% 
      data.frame(check.names = F)
    gmt.path.meta <- read_tsv('https://rest.kegg.jp/list/pathway/map', col_names = F) %>% 
      data.frame(check.names = F)
    cpd.anot <- read_tsv('https://rest.kegg.jp/conv/compound/pubchem', col_names = F) %>% 
      data.frame(check.names = F)
    cpd.anot2 <- read_tsv('https://rest.kegg.jp/conv/compound/chebi', col_names = F) %>% 
      data.frame(check.names = F)
  }, error = function(e) {
    cat("网络不佳，请调整您的网络或使用代理工具\n")
    cat("Error occurred: ", conditionMessage(e), "\n")
  }, warning = function(w) {
    cat("网络不佳，请调整您的网络或使用代理工具\n")
    cat("Warning occurred: ", conditionMessage(w), "\n")
  }, finally = {})
  # data clean 1: original -----------
  names(gmt.path.cpd) <- c("path_map.id", "cpd_kegg.id")
  names(cpd.meta) <- c("cpd_kegg.id", "cpd.name")
  names(gmt.path.meta) <- c("path_map.id", "path_name")
  names(cpd.anot) <- c("cpd_pubchem.id", "cpd_kegg.id")
  names(cpd.anot2) <- c("cpd_chebi.id", "cpd_kegg.id")
  
  # 
  gmt.path.cpd[,1] <- sapply(strsplit(gmt.path.cpd[,1], ":"), function(x) x[[2]])
  gmt.path.cpd[,2] <- sapply(strsplit(gmt.path.cpd[,2], ":"), function(x) x[[2]])
  cpd.anot[,2] <- sapply(strsplit(cpd.anot[,2], ":"), function(x) x[[2]])
  cpd.anot2[,2] <- sapply(strsplit(cpd.anot2[,2], ":"), function(x) x[[2]])
  
  keggpathway <- merge(gmt.path.meta, gmt.path.cpd, by = "path_map.id", all.x = T)
  keggpathway <- merge(keggpathway, cpd.meta, by = "cpd_kegg.id", all.x = T)
  keggpathway <- keggpathway[,c(2,3,1,4)] %>% dplyr::arrange(path_map.id)
  
  # 加上其它id
  dx.df <- merge(keggpathway, cpd.anot, by = "cpd_kegg.id", all.x = T)
  dx.df <- merge(dx.df, cpd.anot2, by = "cpd_kegg.id", all.x = T)
  dx.df <- dx.df %>% 
    dplyr::select(path_map.id, path_name, everything()) %>% 
    dplyr::arrange(path_map.id)
  
  # data clean 2: to user friendly -----------
  allTogmt <- function(dx.df, id = cpd_kegg.id) {
    df <- dx.df %>% dplyr::select(path_map.id, path_name, all_of(id)) %>% 
      rename(term = path_name, gene = id) %>% 
      group_by(term, gene) %>% 
      distinct(.keep_all = T) %>% 
      na.omit()
    return(df)
  }
  df1 <- allTogmt(dx.df = dx.df, id = "cpd_kegg.id")
  df2 <- allTogmt(dx.df = dx.df, id = "cpd_pubchem.id")
  df3 <- allTogmt(dx.df = dx.df, id = "cpd_chebi.id")
  
  metabolismKEGGdb <- list(cpd = df1,
                           chebi = df3,
                           pubmed = df2)
  metabolismKEGGgmt <- lapply(metabolismKEGGdb, function(x) x[,2:3])
  
  original_data <- list(gmt.path.cpd = gmt.path.cpd,
                        cpd.meta = cpd.meta,
                        gmt.path.meta = gmt.path.meta,
                        cpd.anot = cpd.anot,
                        cpd.anot2 = cpd.anot2)
  return(list(keggdf.all = keggpathway, 
              KEGGPATHWAYdb = metabolismKEGGdb,
              KEGGPATHWAYgmt = metabolismKEGGgmt,
              original_data = original_data))
}
