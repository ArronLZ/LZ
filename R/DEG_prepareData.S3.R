#' S3 general function DEG_prepareData()
#' @description S3 general function DEG_prepareData(), prepare the data that LZ diff analysis need,\cr
#' return a list(eset=xx, group=xx, f_mark="xx")\cr
#' DEG_prepareData.default default method \cr
#' DEG_prepareData.DEeset  DEeset method\cr
#' DEG_prepareData.XENA xena method
#'
#' @export
#' @author Jiang
#' @examples
#' # read from XENABD object, oop  ============
#' # DEeset <- DEG_prepareData(XENADB,
#' #                           id_dot = T, col.by = "Ensembl_ID",
#' #                           col.del=NULL, auto.del.character=T,
#' #                           annot_trans=T, f_mark="diff",
#' #                           oop = T,
#' #                           oop.group.suffix1 = "01A",
#' #                           oop.group.suffix2 = "11A",
#' #                           oop.group.endT = "Tumor",
#' #                           oop.group.endF = "Normal")
#' #
#' # read from local file, oop ============
#' # test <- DEG_prepareData(eset_file="D:/TCGA.counts.tsv.gz",
#' #                         id_dot = T, col.by = "Ensembl_ID",
#' #                         col.del=NULL, auto.del.character=T,
#' #                         annot_trans=T, f_mark="diff",
#' #                         oop = T,
#' #                         oop.group.suffix1 = "01A",
#' #                         oop.group.suffix2 = "11A",
#' #                         oop.group.endT = "Tumor",
#' #                         oop.group.endF = "Normal")
#' #
#' # read frome local file(eset, group), no oop ============
#' # 仅当oop = F, 才会进入非oop流程, 此时需要指定group_flie参数才行。默认情况为F
#' # glist <- DEG_prepareData(eset_file = "data/eset.csv", #表达数据的相对路径
#' #                          group_file = "data/group.csv", #分组文件的相对路径
#' #                          id_dot = F,  # ESEM是否有点，有点设为T
#' #                          col.by = "ID",  # 基因名列的列名
#' #                          annot_trans = F, # 是否要注释，如果是EsembleID就需要设置为T
#' #                          f_mark = mark)
#' #' #' # exprset.group = list(eset, group, f_mark)
#' # ** 如果bacth=T，则group中第二列为batch列。
#' # ** 矩阵数据格式(数值型，整型)
#' # gene  | row1 | row2 | row3  | row4
#' # gene1 |  34  |  23  |  56   |  23
#' # gene2 |  35  |  23  |  12   |  23
#' # gene3 |  12  |  78  |  78   |  78
#' # ** 分组数据格式：需要组的行名 包含于 表达谱的列名 rownames(group) %in% colname(eset)
#' # Sample   |   Type    | BATCH
#' # rowname1 |   tumor   |  1
#' # rowname2 |   tumor   |  1
#' # rowname3 |   normal  |  2
#' # rowname4 |   normal  |  2
#' # if is.human = F, please set the orgDb and fromType correctly.
DEG_prepareData <- function(...) {
  UseMethod("DEG_prepareData")
}


#' @param eset_file character, the exps data file name, default: "gene_count.csv"
#' @param eset.islog logical, if the eset have been log2, default is F, meaning the eset data is not log2
#' @param id_dot logical, if the id column is ensemble id with dot, ex.ESEMxxxx.3
#' @param col.by character, delete the duplicate by this column
#' @param col.del character, manually delete some column
#' @param auto.del.character logical, if auto delete the column those class is character
#' @param group_file the group data file name, such as: "group.csv", if the param group_file is missing and oop = T this function will use DEeset$updateGroup() method
#' @param annot_trans logical, do the eset need be annot to symbol? if the datafame's first coloumn is not offical symbol, you can set to TRUE.
#' @param f_mark file mark(optional), default: diff
#' @param is.human logical, default is True, if it is True, the param orgDB,fromTyp and toType will be ignore
#' @param orgDb character, annotation db, ex. "org.Mm.eg.db"
#' @param fromType character, eset id column type
#' @param toType character, the annot result
#' @param oop logical, if use oop flow
#' @param oop.group.suffix1 character, default 01A, eset's name end with 01A will use as group Tumor, DEeset$updateGroup() 
#' @param oop.group.suffix2 character, default 01A, eset's name end with 11A will use as group Normal, DEeset$updateGroup() 
#' @param oop.group.group1 character, default Tumor, eset's name end with 01A will use as group Tumor, DEeset$updateGroup() 
#' @param oop.group.group2 character, default Normal, eset's name end with 01A will use as group Normal, DEeset$updateGroup() 
#' @param oop.group.method  number, default is 1, \cr
#' if method is not 1 will use the existed self$group to make eset2
#' 
#' @return lists
#' @export
#' @rdname DEG_prepareData
#' @importFrom data.table fread
#' @importFrom clusterProfiler bitr
#'
#' @author Jiang
DEG_prepareData.default <- function(eset_file="gene_count.csv",
                                    eset.islog = F,
                                    id_dot = F, col.by = "ID",
                                    col.del=NULL, auto.del.character=T,
                                    group_file,
                                    annot_trans=T, f_mark="diff",
                                    is.human = T, orgDb = "org.Mm.eg.db",
                                    fromType = "ENSEMBL", 
                                    toType = c("SYMBOL", "UNIPROT"),
                                    oop = FALSE,
                                    oop.group.suffix1 = "01A",
                                    oop.group.suffix2 = "11A",
                                    oop.group.group1 = "Tumor",
                                    oop.group.group2 = "Normal",
                                    oop.group.method = 1
) {
  eset <- data.table::fread(eset_file, data.table = F)
  if (eset.islog) { eset[, 2:ncol(eset)] <- round(2^eset[, 2:ncol(eset)] - 1) }
  eset <- quchong(eset = eset, col.by = col.by, col.del = col.del, 
                  auto.del.character = auto.del.character)
  
  # if the id contain ., only keep the part before .  ---------
  if (id_dot) {
    eset$annot <- sapply(strsplit(rownames(eset), "\\."), 
                         function(x) x[[1]])
    eset <- quchong(eset = eset, col.by = "annot")
  }
  # if the id contain ., only keep the part before .///  ---------
  
  # id annot or not  ----------
  if (is.human) {
    # human
    tmpenv <- new.env()
    data(gencode.v22.annot, package = "LZ", envir = tmpenv)
    all_anot <- tmpenv$gencode.v22.annot
    if (annot_trans) {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$ensembl_gene_id, 
      ]
      eset <- cbind(SYMBOL = all_anot$hugo_mRNA[match(rownames(eset), 
                                                      all_anot$hugo_mRNA$ensembl_gene_id), 3], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    } else {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$symbol, 
      ]
    }
    rm(tmpenv)
  } else {
    # not human , default is mm
    if (annot_trans) {
      annot.df.org <- bitr(rownames(eset), fromType = fromType, 
                           toType = toType, OrgDb = orgDb)
      eset <- eset[rownames(eset) %in% annot.df.org[, 1], ]
      eset <- cbind(SYMBOL = annot.df.org[match(rownames(eset), 
                                                annot.df.org[, 1]), "SYMBOL"], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    }
  }
  
  
  # sekf function
  checkgroup <- function(eset, group) {
    # check if the rowname of group is all from the colnames of eset -----
    if ( all(rownames(group) %in% names(eset)) ) {
      # eset's coloumns > group's rows
      if (length(names(eset)) > length(rownames(group))) {
        cat(" 请注意：!!!\n 提供的分组文件中样本数量少于表达矩阵中的样本数：\n
 这意味着将会只保留分组文件中的样本来进行后续分析，\n 不在分组文件中的样本将会被剔除!!!\n
 分组文件内容展示：\n")
        print(group)
      }
      # only keep the eset's coloumn by group's rownames sample
      eset <- eset[, rownames(group)]
      # check all eset's coloumn == group's rowname
      if (all(names(eset) == rownames(group))) {
        glist <- list(eset = eset, group = group, f_mark = f_mark)
        glist$group[,1] <- factor(glist$group[,1])
      } else {
        stop("提供的分组信息和表达信息不匹配，请检查后再运行\n")
      }
    } else {
      stop("提供的分组信息和表达信息不匹配(分组文件有样本名不在表达矩阵的样本名中)，\n
         请检查后再运行\n")
    }
    return(glist)
  }
  
  
  if (!oop) {
    # id annot or not ///  ----------
    group <- read.csv(group_file, row.names = 1)
    glist <- checkgroup(eset = eset, group = group)
    # check if the rowname of group is all from the colnames of eset /// -----
    return(glist)
  } 
  
  # OOP ---------
  if (oop) {
    if (missing(group_file)) {
      group <- NULL
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = eset)
      DEeset$updateGroup_Eset(
        suffix1 = oop.group.suffix1,
        suffix2 = oop.group.suffix2,
        group1 = oop.group.group1,
        group2 = oop.group.group2,
        method = oop.group.method
      )
      return(DEeset)
    } else {
      group <- read.csv(group_file, row.names = 1)
      glist <- checkgroup(eset = eset, group = group)
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = glist$eset, group = glist$group)
      return(DEeset)
    }
  }
}

#' @param obj obeject the object DEeset
#' @param f_mark character a analysis mark(optional), default: "DE"
#'
#' @return list
#' @export
#' @rdname DEG_prepareData
#' @method DEG_prepareData DEeset
DEG_prepareData.DEeset <- function(obj, f_mark = "DE") {
  glist <- list(eset = obj$eset, group = obj$group, f_mark = f_mark)
  return(glist)
}


#' @return list
#' @export
#' @rdname DEG_prepareData
DEG_prepareData.XENA <- function(obj, # eset_file="gene_count.csv"
                                 eset.islog = T,
                                 id_dot = T, col.by = "Ensembl_ID",
                                 col.del=NULL, auto.del.character=T,
                                 group_file, 
                                 annot_trans=T, f_mark="diff",
                                 is.human = T, orgDb = "org.Mm.eg.db",
                                 fromType = "ENSEMBL", 
                                 toType = c("SYMBOL", "UNIPROT"),
                                 oop = T,
                                 oop.group.suffix1 = "01A",
                                 oop.group.suffix2 = "11A",
                                 oop.group.group1 = "Tumor",
                                 oop.group.group2 = "Normal",
                                 oop.group.method = 1
                                 ) {
  eset <- obj$eset.count
  if (eset.islog) { eset[, 2:ncol(eset)] <- round(2^eset[, 2:ncol(eset)] - 1) }
  eset <- quchong(eset = eset, col.by = col.by, col.del = col.del, 
                  auto.del.character = auto.del.character)
  
  # if the id contain ., only keep the part before .  ---------
  if (id_dot) {
    eset$annot <- sapply(strsplit(rownames(eset), "\\."), 
                         function(x) x[[1]])
    eset <- quchong(eset = eset, col.by = "annot")
  }
  # if the id contain ., only keep the part before .///  ---------
  
  # id annot or not  ----------
  if (is.human) {
    # human
    tmpenv <- new.env()
    data(gencode.v22.annot, package = "LZ", envir = tmpenv)
    all_anot <- tmpenv$gencode.v22.annot
    if (annot_trans) {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$ensembl_gene_id, 
      ]
      eset <- cbind(SYMBOL = all_anot$hugo_mRNA[match(rownames(eset), 
                                                      all_anot$hugo_mRNA$ensembl_gene_id), 3], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    } else {
      eset <- eset[rownames(eset) %in% all_anot$hugo_mRNA$symbol, 
      ]
    }
    rm(tmpenv)
  } else {
    # not human , default is mm
    if (annot_trans) {
      annot.df.org <- bitr(rownames(eset), fromType = fromType, 
                           toType = toType, OrgDb = orgDb)
      eset <- eset[rownames(eset) %in% annot.df.org[, 1], ]
      eset <- cbind(SYMBOL = annot.df.org[match(rownames(eset), 
                                                annot.df.org[, 1]), "SYMBOL"], eset)
      eset <- quchong(eset, col.by = "SYMBOL")
    }
  }
  
  
  # sekf function
  checkgroup <- function(eset, group) {
    # check if the rowname of group is all from the colnames of eset -----
    if ( all(rownames(group) %in% names(eset)) ) {
      # eset's coloumns > group's rows
      if (length(names(eset)) > length(rownames(group))) {
        cat(" 请注意：!!!\n 提供的分组文件中样本数量少于表达矩阵中的样本数：\n
 这意味着将会只保留分组文件中的样本来进行后续分析，\n 不在分组文件中的样本将会被剔除!!!\n
 分组文件内容展示：\n")
        print(group)
      }
      # only keep the eset's coloumn by group's rownames sample
      eset <- eset[, rownames(group)]
      # check all eset's coloumn == group's rowname
      if (all(names(eset) == rownames(group))) {
        glist <- list(eset = eset, group = group, f_mark = f_mark)
        glist$group[,1] <- factor(glist$group[,1])
      } else {
        stop("提供的分组信息和表达信息不匹配，请检查后再运行\n")
      }
    } else {
      stop("提供的分组信息和表达信息不匹配(分组文件有样本名不在表达矩阵的样本名中)，\n
         请检查后再运行\n")
    }
    return(glist)
  }
  
  
  if (!oop) {
    # id annot or not ///  ----------
    group <- read.csv(group_file, row.names = 1)
    glist <- checkgroup(eset = eset, group = group)
    # check if the rowname of group is all from the colnames of eset /// -----
    return(glist)
  } 
  
  # OOP ---------
  if (oop) {
    if (missing(group_file)) {
      group <- NULL
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = eset)
      DEeset$updateGroup_Eset(
        suffix1 = oop.group.suffix1,
        suffix2 = oop.group.suffix2,
        group1 = oop.group.group1,
        group2 = oop.group.group2,
        method = oop.group.method
      )
      return(DEeset)
    } else {
      group <- read.csv(group_file, row.names = 1)
      glist <- checkgroup(eset = eset, group = group)
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = glist$eset, group = glist$group)
      return(DEeset)
    }
  }
}


