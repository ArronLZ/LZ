#' Prepare data for DEG anlysis
#' @description S3 general function DEG_prepareData(), prepare the data that LZ diff analysis need,\cr
#' return a list(eset=xx, group=xx, f_mark="xx")
#'
#' @param object any R object
#' @param eset_file character, the exps data file name, default: "gene_count.csv"
#' @param eset.islog logical, if the eset have been log2, default is F, meaning the eset data is not log2
#' @param id_dot logical, if the id column is ensemble id with dot, ex.ESEMxxxx.3
#' @param col.by character, delete the duplicate by this column
#' @param col.del character, manually delete some column
#' @param auto.del.character logical, if auto delete the column those class is character
#' @param annot_trans logical, do the eset need be annot to symbol? if the datafame's first coloumn is not offical symbol, you can set to TRUE.
#' @param group_file the group data file name, such as: "group.csv", if the param group_file is missing and oop = T this function will use DEeset$updateGroup() method
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
#' @param ... other param
#'
#' @export
#' @importFrom data.table fread
#' @importFrom clusterProfiler bitr
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
DEG_prepareData <- function(object, eset_file, eset.islog, id_dot, col.by,
                            col.del=NULL, auto.del.character=T, annot_trans,
                            group_file, f_mark="diff",
                            is.human = T,
                            orgDb = "org.Mm.eg.db", fromType = "ENSEMBL",
                            toType = c("SYMBOL", "UNIPROT"),
                            oop = FALSE,
                            oop.group.suffix1, oop.group.suffix2,
                            oop.group.group1, oop.group.group2, ...) {
  UseMethod("DEG_prepareData")
}


#' @method DEG_prepareData default
#' @rdname DEG_prepareData
#' @export
DEG_prepareData.default <- function(...,
                                    eset_file="gene_count.csv",
                                    eset.islog = F,
                                    id_dot = F, col.by = "ID",
                                    col.del=NULL, auto.del.character=T,
                                    annot_trans=T,
                                    group_file,
                                    f_mark="diff",
                                    is.human = T, orgDb = "org.Mm.eg.db",
                                    fromType = "ENSEMBL",
                                    toType = c("SYMBOL", "UNIPROT"),
                                    oop = FALSE,
                                    oop.group.suffix1 = "01A",
                                    oop.group.suffix2 = "11A",
                                    oop.group.group1 = "Tumor",
                                    oop.group.group2 = "Normal") {
  stopifnot(!missing(eset_file))
  if (!oop) {
    if ( missing(group_file) ) {
      stop("如果oop = F，则必须设置group_file参数")
    }
  }
  
  # eset <- data.table::fread(eset_file, data.table = F)
  if (is.data.frame(eset_file) | is.character(eset_file)) {
    tryCatch({
      if (is.data.frame(eset_file)) {
        eset <- eset_file
      } else {
        eset <- data.table::fread(eset_file, data.table = F)
      }
    }, error = function(e) {
      stop("eset_file must be data.frame or character(a valid file path)")
    })
  } else {
    stop("eset_file must be data.frame or character")
  }
  
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
  
  
  # sekf function 1
  read_group <- function() {
    # group <- read.csv(group_file, row.names = 1)
    if (is.data.frame(group_file) | is.character(group_file)) {
      tryCatch({
        if (is.data.frame(group_file)) {
          group <- group_file
        } else {
          group <- read.csv(group_file, row.names = 1)
        }
      }, error = function(e) {
        stop("eset_file must be data.frame or character(a valid file path)")
      })
    } else {
      stop("eset_file must be data.frame or character")
    }
    return(group)
  }
  # sekf function 2
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
    group <- read_group()
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
        method = 1)
      return(DEeset)
    } else {
      group <- read_group()
      glist <- checkgroup(eset = eset, group = group)
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = glist$eset, group = glist$group)
      DEeset$updateGroup_Eset(
        #suffix1 = oop.group.suffix1,
        #suffix2 = oop.group.suffix2,
        #group1 = oop.group.group1,
        #group2 = oop.group.group2,
        method = 2)
      return(DEeset)
    }
  }
}


#' @method DEG_prepareData XENA
#' @rdname DEG_prepareData
#' @export
DEG_prepareData.XENA <- function(object, ..., eset.islog = T, id_dot = T,
                                 col.by = "Ensembl_ID", col.del=NULL,
                                 auto.del.character=T, annot_trans=T,
                                 group_file,
                                 f_mark="diff",
                                 is.human = T,
                                 orgDb = "org.Mm.eg.db",
                                 fromType = "ENSEMBL",
                                 toType = c("SYMBOL", "UNIPROT"),
                                 oop = T,
                                 oop.group.suffix1 = "01A",
                                 oop.group.suffix2 = "11A",
                                 oop.group.group1 = "Tumor",
                                 oop.group.group2 = "Normal"
) {
  if (!oop) {
    if ( missing(group_file) ) {
      stop("如果oop = F，则必须设置group_file参数")
    }
  }
  eset <- object$eset.count
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
        method = 1)
      return(DEeset)
    } else {
      group <- read.csv(group_file, row.names = 1)
      glist <- checkgroup(eset = eset, group = group)
      DEeset <- YZ::DEeset$new(mark = f_mark, eset = glist$eset, group = glist$group)
      DEeset$updateGroup_Eset(
        suffix1 = oop.group.suffix1,
        suffix2 = oop.group.suffix2,
        group1 = oop.group.group1,
        group2 = oop.group.group2,
        method = 2)
      return(DEeset)
    }
  }
}


#' @method DEG_prepareData DEeset
#' @rdname DEG_prepareData
#' @export
DEG_prepareData.DEeset <- function(object, ..., f_mark = "DE") {
  glist <- list(eset = object$eset2, group = object$group, f_mark = f_mark)
  return(glist)
}


#' Create a new XENA object
#' @description Create a new XENA object
#' 
#' @param eset.count dataframe, the exprs dataframe, column is gene_id\cr
#' the last column is the sample, row is the gene.
#' @param info character, the information of eset.count
#'
#' @return a XENA object
#' @export
newXENA <- function(eset.count, info) {
  # requireNamespace("YZ")
  xenadb <- YZ::Class("XENA") #YZ package
  xenadb$new(xenadb, "eset.count", eset.count, info = info)
  return(xenadb)
}


#' Create a new DEeset object
#' @description Create a new DEeset object by exsited eset and group object
#' 
#' @param eset dataframe, the cleaned exprs dataframe, rowname is gene_id
#' @param group dataframe, the cleaned group dataframe, rowname is sample, \cr
#' the first cloumn is the Type and must be factor
#' @param f_mark character, the mark of this DEeset object, default is 'DE'
#'
#' @return a DEeset object
#' @export
#'
#' @author Jiang
newDEeset <- function(eset, group, f_mark = "DE") {
  # requireNamespace("YZ")
  stopifnot(all(colnames(eset) == rownames(group)),
            is.factor(group[,1]),
            names(group)[1] == "Type"
            )
  # group <- read.csv(group_file, row.names = 1)
  # glist <- checkgroup(eset = eset, group = group)
  DEeset <- YZ::DEeset$new(mark = f_mark, 
                           eset = eset, 
                           group = group)
  DEeset$eset2 <- eset
  return(DEeset)
}
