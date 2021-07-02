## RNAseq 工作流程

    library(RNAseqFlow)
    input = createCountPhe()

    dds = DEseqObj(count_data=input[[1]],col_data=input[[2]],design_names = "condition+type",group_name ="condition",ref_level = "untreated")
    #> [1] "The id order between gene count file and phenotype file is identical without modification!"
    res = DESeqRes(dds,foldChange=0.58,adjPvalue=0.05)
    DESeqObjPCA(dds) # plot PCA

![](../Figs/unnamed-chunk-2-1.png)

    DESeqResVolcano(res) # plot volcano

![](../Figs/unnamed-chunk-2-2.png)

    DEsig = row.names(res)[res$regulate!="Normal"]

## 富集分析

加载一些程序包：

    library(org.Dm.eg.db)
    library(clusterProfiler)

### GO 过度代表检验 (over-representation test)

    goRes = GO(DEsig,type="ENSEMBL",db=org.Dm.eg.db)

为了方便查看，需要把GO和KEGG的结果保存为文件时，base\_name参数可以实现这个功能。

### KEGG 过度代表检验 (over-representation test)

    keggRes = KEGG(DEsig,type="ENSEMBL",organism="dme",db=org.Dm.eg.db);
    head(keggRes)

值得注意的是，当进行KEGG或gseKEGG时，如果物种是人类，那么基因名为ENTREZID时，对应的keyType是“kegg”或“ncbi-geneid”。当对其它物种进行分析时，ENTREZID对应的是ncbi-geneid。

### GO的GSEA分析

    geneList = res$log2FoldChange
    names(geneList)=row.names(res)
    gseagoRes = GSEAgo(geneList,type="ENSEMBL",db=org.Dm.eg.db)
    head(gseagoRes)
    if(nrow(gseagoRes)>=1){
      print(gseaplot(gseagoRes,geneSetID=gseagoRes$ID[1],title=paste("BP : ",gseagoRes$Description[1],sep="")))
    }

### KEGG的GSEA分析

    gseakeggRes = GSEAkegg(geneList,type="ENSEMBL",organism = "dme",db=org.Dm.eg.db,pvalueCutoff=1)


    if(nrow(gseakeggRes)>1){
      print(gseaplot(gseakeggRes,geneSetID=gseakeggRes$ID[1],title=paste("KEGG : ",gseakeggRes$Description[1],sep="")))
    }

![](../Figs/unnamed-chunk-7-1.png)

## Todo

GSEAGO出错
