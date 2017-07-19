context("conversion from DNAbin to genind")


test_that("genind objects can be converted to DNAbin", {
  skip_if_not_installed("ape")
  data("woodmouse", package = "ape")
  wm <- as.genind.DNAbin(woodmouse, rep(c("A", "B", "C"), each=5))
  expect_is(wm, "genind")
  expect_equal(nrow(woodmouse), nInd(wm))
})
