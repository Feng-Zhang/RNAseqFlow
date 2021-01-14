##' @title GO
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param ids a character vector of gene id
##' @param type the type of ids which must be one of idType(OrgDb = db) shows
##' @param db the database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
##' @param outputDir the path of output
##' @param baseName the prefix name of outputs would be saved
##' @return NULL
##' @importFrom stringr str_split_fixed
##' @importFrom clusterProfiler bitr
##' @importFrom utils write.csv
##' @importFrom grDevices pdf png dev.off
##' @importFrom graphics barplot
##' @export
##'

GO = function(ids,type="SYMBOL",db="org.Hs.eg.db",outputDir=".",baseName="DEG"){
  if(type!="ENTREZID") ids = bitr(ids, fromType=type, toType="ENTREZID", OrgDb=db)[,2]

  ego <- enrichGO(gene= ids,#括号内必须加逗号，否则是矩阵
                  OrgDb= db,
                  ont= "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,  readable = TRUE)
  write.csv(ego,file=paste0(outputDir,"/",baseName,"_GO_enrichment.csv"),row.names = F)


  print("############# barplot and dot plot for go enrichment ####################")
  pdf(file = paste0(outputDir,"/",baseName,"_GO_barplot.pdf"),height=4,width=8)
  if(!is.null(ego)) print(barplot(ego, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()
  png(filename = paste0(outputDir,"/",baseName,"_GO_barplot.png"),height=2000,width=5000,res=300,type="cairo")
  if(!is.null(ego)) print(barplot(ego, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()


  pdf(file = paste0(outputDir,"/",baseName,"_GO_dotplot.pdf"),height=4,width=8)
  if(!is.null(ego)) print(dotplot(ego,showCategory=20)) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()

  png(filename = paste0(outputDir,"/",baseName,"_GO_dotplot.png"),height=2000,width=5000,res=300,type="cairo")
  if(!is.null(ego)) print(dotplot(ego,showCategory=20)) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()
}



