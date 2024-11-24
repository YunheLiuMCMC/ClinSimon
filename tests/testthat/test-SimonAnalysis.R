test_that("Test Post-Trial Inference work correctly", {
  # Case 1 for ATSS Simon Design in using `MIDp` confidence interval option
  Inference.case1 <- c(0.322, 0.108, 0.538, 0.831)
  expect_equal(as.numeric(SimonAnalysis(2,7,13,30,3,12,0.05,c(0.025,0.975),"MIDp",0.40)),
               Inference.case1)

  # Case 2 for ATSS Simon Design in using `Jung` confidence interval option
  Inference.case2 <- c(0.429, 0.257, 0.568, 0.019)
  expect_equal(as.numeric(SimonAnalysis(2,16,11,28,2,13,0.077,c(0.025,0.975),"Jung",0.25)),
               Inference.case2)
})
