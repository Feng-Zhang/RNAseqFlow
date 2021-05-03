##' @title GO
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param ids a character vector of gene id
##' @param type the type of ids which must be one of idType(OrgDb = db) shows
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
##' @param db the annotation database of organism, such as org.Hs.eg.db for human
##' @param outputDir the path of output
##' @param baseName the prefix name of outputs would be saved
##' @return NULL
##' @importFrom stringr str_split_fixed
##' @importFrom clusterProfiler bitr enrichKEGG
##' @importFrom utils write.csv
##' @importFrom grDevices pdf png dev.off
##' @importFrom graphics barplot
##' @importFrom enrichplot dotplot
##' @importFrom org.Hs.eg.db org.Hs.eg.db
##' @export
##'

KEGG = function(ids,type="SYMBOL",organism="hsa",db=org.Hs.eg.db,outputDir=".",baseName="Test"){
  dgeKegg <- NULL
  if(type!="ENTREZID") ids = bitr(ids, fromType=type, toType="ENTREZID", OrgDb=db)[,2]
  print("#####################  barplot and dot plot for KEGG pathway ################################")
  dgeKegg = enrichKEGG(gene = ids, organism=organism, pvalueCutoff = 0.05)
  write.csv(dgeKegg,file=paste0(outputDir,"/",baseName,"_KEGG_enrichment.csv"),row.names = F)

  pdf(file = paste0(outputDir,"/",baseName,"_KEGG_dotplot.pdf"),height=4,width=8)
  if(!is.null(dgeKegg)) print(dotplot(dgeKegg,showCategory=20))  else plot(0,main="Sorry,No gene can be enriched")
  dev.off()

  png(filename = paste0(outputDir,"/",baseName,"_KEGG_dotplot.png"),height=2000,width=5000,res=300,type="cairo")
  if(!is.null(dgeKegg)) print(dotplot(dgeKegg,showCategory=20))  else plot(0,main="Sorry,No gene can be enriched")
  dev.off()

  pdf(file = paste0(outputDir,"/",baseName,"_KEGG_barplot.pdf"),height=4,width=8)
  if(!is.null(dgeKegg)) print(barplot(dgeKegg, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()

  png(filename = paste0(outputDir,"/",baseName,"_KEGG_barplot.png"),height=2000,width=5000,res=300,type="cairo")
  if(!is.null(dgeKegg)) print(barplot(dgeKegg, showCategory=20, x = "GeneRatio")) else plot(0,main="Sorry,No gene can be enriched")
  dev.off()
}
