##' @title GSEAgo
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param geneChangeList a numeric vector of the fold change of genes. It must have name with ENTREZID.
##' @param db the database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
##' @return NULL
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' ego = GSEAgo(geneList)
##' }
##' @importFrom clusterProfiler gseGO  setReadable
##' @export
##'

GSEAgo = function(geneChangeList,db=org.Hs.eg.db){
  #geneChangeList=geneList
  ego <- NULL
  geneChangeList = sort(geneChangeList, decreasing = TRUE)
  ego = gseGO(geneChangeList, ont="BP",OrgDb=db,eps=0)
  ego = setReadable(ego, OrgDb =db, keyType="ENTREZID")
  # if(nrow(gseaGo)>1){
  #   write.csv(gseaGo,file=paste(trait,"_gseaGo.csv",sep=""))
  #   png(paste(trait,"_gseaGo",".png",sep=""),height=1000,width=1500,res=150,type="cairo")
  #   print(gseaplot(gseaGo,geneSetID=gseaGo$ID[1],title=paste("BP : ",gseaGo$Description[1],sep="")))
  #   dev.off()
  # }
  return(ego)
}


##' @title GSEAkegg
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param geneChangeList a numeric vector of the fold change of genes. It must have name with ENTREZID.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html', such as hsa for human, mmu for mouse.
##' @param db the annotation database of specific organism
##' @return NULL
##' @importFrom clusterProfiler gseKEGG setReadable
##' @importFrom org.Hs.eg.db org.Hs.eg.db
##' @export
##'

GSEAkegg = function(geneChangeList,organism="hsa",db=org.Hs.eg.db){
  ekegg <- NULL
  ## gsea
  geneChangeList = sort(geneChangeList, decreasing = TRUE)
  ekegg = gseKEGG(geneChangeList,organism=organism,pvalueCutoff = 0.05,eps=0)
  ekegg = setReadable(ekegg, OrgDb = db, keyType="ENTREZID")
  #
  # if(nrow(gseaKegg)>1){
  # write.csv(gseaKegg,file=paste(trait,"_gseaKegg.csv",sep=""))
  # png(paste(trait,"_gseaKegg",".png",sep=""),height=1000,width=1500,res=150,type="cairo")
  # print(gseaplot(gseaKegg,geneSetID=gseaKegg$ID[1],title=paste("KEGG : ",gseaKegg$Description[1],sep="")))
  # dev.off()
  # }
  return(ekegg)
}

