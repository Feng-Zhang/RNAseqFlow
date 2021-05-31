##' @title DESeqObjPCA
##' @description This function would plot PCA based on DESeqDataSet object and save as pdf and png format
##' @param dds a matrix with row name
##' @param group_name the number of colData to split case and control group
##' @param outputDir the path of output
##' @param base_name NULL or character, the prefix of output. When value is NULL, the figures are not saved.
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



DESeqObjPCA = function(dds,group_name="condition",outputDir=".",base_name=NULL){
  #groupNum=1;
  name <- NULL
  stopifnot("DESeqDataSet"%in% class(dds))
  if(nrow(dds)>1000) vsd=vst(dds, blind=FALSE) else vsd=varianceStabilizingTransformation(dds)
  p=plotPCA(vsd, intgroup=group_name)+geom_text(aes(label=name), vjust = 'inward', hjust = 'inward')

  if(!is.null(base_name)){
    jpeg(filename=paste0(outputDir,"/",base_name,"_PCA.jpg"),res=300,height = 3000,width=3000,type="cairo")
    print(p)
    dev.off()

    pdf(file=paste0(outputDir,"/",base_name,"_PCA.pdf"))
    print(p)
    dev.off()
  }
  return(p)
}
