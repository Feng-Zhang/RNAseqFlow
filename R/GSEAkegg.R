##' @title GSEAkegg
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param geneChangeList a numeric vector of the fold change of genes. It must have name with ENTREZID.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
##' @param outputDir the path of output
##' @param baseName the prefix name of outputs would be saved
##' @return NULL
##' @importFrom clusterProfiler gseKEGG setReadable
##' @export
##'

GSEAkegg = function(geneChangeList,organism="hsa",outputDir=".",baseName="DEG"){

  ## gsea
  gseaKegg = gseKEGG(geneChangeList,organism="hsa",pvalueCutoff = 0.05,eps=0)
  gseaKegg = setReadable(gseaKegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  #
  # if(nrow(gseaKegg)>1){
  # write.csv(gseaKegg,file=paste(trait,"_gseaKegg.csv",sep=""))
  # png(paste(trait,"_gseaKegg",".png",sep=""),height=1000,width=1500,res=150,type="cairo")
  # print(gseaplot(gseaKegg,geneSetID=gseaKegg$ID[1],title=paste("KEGG : ",gseaKegg$Description[1],sep="")))
  # dev.off()
  # }
  return(gseakegg)
}
