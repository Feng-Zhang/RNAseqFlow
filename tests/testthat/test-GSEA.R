context("Test enrichment analysis including GO, KEGG and GSEA")
#library(DOSE)

test_that("Test GSEA for GO and KEGG term", {
  data(geneList,package = "DOSE")
  expect_error(GSEAgo(geneList,db = "org.Hs.eg.db"),"unable to find an inherited method")
  ego = GSEAgo(geneList)
  expect_is(ego,"gseaResult")
  # egoKegg = GSEAkegg(geneList)
})


