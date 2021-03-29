context("Create test data including count and phenotype data")

test_that("input is a list", {
  input = createCountPhe()
  countData = input[[1]]
  colData = input[[2]]
  expect_true(any(class(countData)%in%c("matrix","array")))
  expect_equal(class(colData),"data.frame")
})
