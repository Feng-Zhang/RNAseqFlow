##' @title GSEAgo
##' @description This function would plot volcano based on DESeqDEGres output and save as pdf and png format
##' @param geneChangeList a numeric vector of the fold change of genes. It must have name with ENTREZID.
##' @param db the database, such as org.Hs.eg.db , org.Mm.eg.db and so on.
##' @param outputDir the path of output
##' @param baseName the prefix name of outputs would be saved
##' @return NULL
##' @importFrom clusterProfiler gseGO  setReadable
##' @export
##'

GSEAgo = function(geneChangeList,db="org.Hs.eg.db",outputDir=".",baseName="DEG"){

  geneChangeList = sort(geneChangeList, decreasing = TRUE)
  gseaGo = gseGO(geneChangeList, ont="BP",OrgDb=db,eps=0)
  gseaGo = setReadable(gseaGo, OrgDb =db, keyType="ENTREZID")

  # if(nrow(gseaGo)>1){
  #   write.csv(gseaGo,file=paste(trait,"_gseaGo.csv",sep=""))
  #   png(paste(trait,"_gseaGo",".png",sep=""),height=1000,width=1500,res=150,type="cairo")
  #   print(gseaplot(gseaGo,geneSetID=gseaGo$ID[1],title=paste("BP : ",gseaGo$Description[1],sep="")))
  #   dev.off()
  # }
  return(gseaGO)
}

