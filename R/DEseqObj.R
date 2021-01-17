##' @title DEseqObj
##' @description This function would generate a object of DESeq2 class
##' @param countData a matrix with row name
##' @param colData a dataframe include group information
##' @param groupNum the number of colData to split case and control group
##' @param refLevel a character which is the reference group
##' @return a list
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],refLevel="untreated")
##' }
##' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
##' @importFrom stats relevel as.formula
##' @export
##'
DEseqObj = function(countData,colData,groupNum=1,refLevel="control"){
  #countData=input[[1]];colData=input[[2]];groupNum=1
  if(all(rownames(colData) == colnames(countData))){
    print("The id order between gene count file and phenotype file is identical without modification!")
  } else {
    countData <- countData[, rownames(colData)]
    print(paste0("After modifying, the id order between gene count file and phenotype file is ",all(rownames(colData) == colnames(countData))))
  }

  # analysis and plot
  groupName=colnames(colData)[groupNum]
  # gene expression analysis
  testFormula=paste0("~ ",groupName)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(testFormula) )
  try(`$`(dds,groupName) <- relevel(`$`(dds,groupName), ref = refLevel),silent=TRUE)
  dds <- dds[ rowSums(counts(dds)) > 1 & apply(counts(dds),1,median)>0, ]
  dds <- DESeq(dds)

  return(dds)
}
