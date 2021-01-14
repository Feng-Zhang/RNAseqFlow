##' @title createCountPhe
##' @description This function would generate a list including a countData and colData as input data for testing
##' @return a list
##' @examples input = createCountPhe()
##' @import pasilla
##' @importFrom utils read.csv
##' @importFrom stats median
##' @export

createCountPhe = function(){
  #library("pasilla")

  # import expression data
  pasCts <- system.file("extdata",
                        "pasilla_gene_counts.tsv",
                        package="pasilla", mustWork=TRUE)
  countData <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))

  # import phenotype data
  pasAnno <- system.file("extdata",
                         "pasilla_sample_annotation.csv",
                         package="pasilla", mustWork=TRUE)
  colData <- read.csv(pasAnno, row.names=1)
  colData <- colData[,c("condition","type")]
  colData$condition <- factor(colData$condition)
  colData$type <- factor(colData$type)
  rownames(colData) <- sub("fb", "", rownames(colData))

  # check if the order of individuals between expression and phenotype data
  all(rownames(colData) %in% colnames(countData))
  all(rownames(colData) == colnames(countData))
  countData <- countData[, rownames(colData)]
  all(rownames(colData) == colnames(countData))

  # clean expression data by keeping genes with mean count >0, median count >0
  countData = countData[rowMeans(countData) > 1 & apply(countData,1,median)>0 ,]

  return(list(countData,colData))
}
