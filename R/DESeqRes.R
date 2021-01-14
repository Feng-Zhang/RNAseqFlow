##' @title DESeqRes
##' @description This function would get the results based on DESeqDataSet object
##' @param dds a matrix with row name
##' @param foldChange a dataframe include group information
##' @param adjPvalue the number of colData to split case and control group
##' @return a dataframe
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],refLevel="untreated")
##' res = DESeqRes(dds)
##' }
##' @importFrom  DESeq2 results
##' @export
##'


DESeqRes = function(dds,foldChange=1,adjPvalue=0.05){
  res <- results(dds)
  res <- res[order(res$padj), ]
  res$regulate <- "Normal"
  res[res$log2FoldChange< -foldChange & res$padj< adjPvalue & !is.na(res$padj), "regulate"] <- "Down"
  res[res$log2FoldChange> foldChange & res$padj < adjPvalue & !is.na(res$padj), "regulate"] <- "Up"
  return(as.data.frame(res))
}
