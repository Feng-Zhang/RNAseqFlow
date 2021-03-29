context("RNA-seq analysis basic workfolw, including DE analysis")
test_that("RNA-seq analysis basic workfolw", {
  input = createCountPhe()
  dds = DEseqObj(input[[1]],input[[2]],refLevel="untreated")
  expect_is(dds,"DESeqDataSet")
  res = DESeqRes(dds)
  expect_is(res,"data.frame")
})
