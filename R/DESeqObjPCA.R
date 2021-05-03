##' @title DESeqObjPCA
##' @description This function would plot PCA based on DESeqDataSet object and save as pdf and png format
##' @param dds a matrix with row name
##' @param groupNum the number of colData to split case and control group
##' @param outputDir the path of output
##' @return NULL
##' @examples
##' \dontrun{
##' input = createCountPhe()
##' dds = DEseqObj(input[[1]],input[[2]],refLevel="untreated")
##' DESeqObjPCA(dds)
##' }
##' @importFrom ggplot2 geom_text
##' @importFrom grDevices pdf png dev.off
##' @importFrom DESeq2 plotPCA vst varianceStabilizingTransformation
##' @export
##'



DESeqObjPCA = function(dds,groupNum=1,outputDir="."){
  #groupNum=1;
  name <- NULL
  if(nrow(dds)>1000) vsd=vst(dds, blind=FALSE) else vsd=varianceStabilizingTransformation(dds)
  groupName = colnames(dds@colData)[groupNum]
  p=plotPCA(vsd, intgroup=groupName)+geom_text(aes(label=name), vjust = 'inward', hjust = 'inward')

  png(filename=paste0(outputDir,"/",groupName,"_PCA.png"),res=300,height = 3000,width=3000,type="cairo")
  print(p)
  dev.off()

  pdf(file=paste0(outputDir,"/",groupName,"_PCA.pdf"))
  print(p)
  dev.off()
}
