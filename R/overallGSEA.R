##' @title overallGSEA
##' @description This function would generate bar plot for GSEA results
##' @param data Default dataset to use for plot.
##' @param x Default x of aesthetic mappings to use for plot. Generally, it is Description for enrichment pathway.
##' @param y Default y of aesthetic mappings to use for plot. Generally, it is NES.
##' @return a ggplot2 object
##' @examples
##' \dontrun{
##'  data(geneList, package="DOSE")
##'  geneChangeList = sort(geneList, decreasing = TRUE)
##'  ego = gseGO(geneChangeList, ont="BP",OrgDb=org.Hs.eg.db,eps=0)
##'  ego = setReadable(ego, OrgDb =org.Hs.eg.db)
##'  data = ego[127:147,1:9]
##'  print(overallGSEA(data,Description,NES))
##'  }
##' @import ggplot2
##' @export

overallGSEA = function(data,x,y){
  xExpr = substitute(x)
  yExpr = substitute(y)
  data = data[order(eval(yExpr,data)),]
  data[,"regulate"] = ifelse(eval(yExpr,data)>0,"Up","Down")
  #  data[,"regulate"] = as.factor(data[,"regulate"])
  data[,deparse(xExpr)] = factor(eval(xExpr,data),levels = eval(xExpr,data))
  p = ggplot(data,aes_(xExpr,yExpr,fill=substitute(regulate)))+
    geom_bar(stat = "identity")+
    coord_flip()+scale_fill_manual(values =c("blue","red"))+
    xlab("") + ylab("Normalized Enrichment Score (NES)")+labs(title = "Overall GSEA")+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0))
  return(p)
}
