##' @title DESeqRes
##' @description This function would get the results based on DESeqDataSet object
##' @param geneCountFile a matrix with row name
##' @param pheFile a dataframe include group information
##' @param outputDir the number of colData to split case and control group
##' @param groupNum the number of colData to split case and control group
##' @param foldChange the significant threshold for log2 fold change
##' @param adjPvalue the significant threshold for adjust P-value
##' @param db the database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
##' @return a dataframe
##' @importFrom  stringr str_split_fixed
##' @importFrom utils read.csv
##' @export
##'



RNAseqFlow = function(geneCountFile,pheFile,outputDir=".",groupNum=1,
                      foldChange=1,adjPvalue=0.05,db="org.Hs.eg.db",organism="hsa"){
  # import data
  countData <- as.matrix(read.csv(geneCountFile, row.names="gene_id",check.names=F))
  colData = read.csv(pheFile,header=TRUE,row.names=1)

  # import count data
  dds = DEseqObj(countData,colData,groupNum=groupNum,refLevel="control")
  groupName = colnames(colData)[groupNum]
  res = DESeqRes(dds)

  DESeqObjPCA(dds) # plot PCA
  DESeqResVolcano(res) # plot volcano

  #plot GO and KEGG
  DEsig = row.names(res)[res$regulate!="Normal"]
  ids = unique(str_split_fixed(DEsig,pattern="\\|",n=2)[,2])
  GO(ids,type="SYMBOL",db=db,outputDir=outputDir,baseName=groupName)
  KEGG(ids,type="SYMBOL",organism=organism,outputDir=outputDir,baseName=groupName)

  #PPI

  # GSEA

  save.image(paste0(outputDir,"/",groupName,".RData")) # save image
}
