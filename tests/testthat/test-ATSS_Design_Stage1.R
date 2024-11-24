test_that("Test ATSS_Design_Stage1() get correct design output", {
  # Case 1 for Adaptive Threshold and Sample Size Simon Design (ATSS Simon) in Stage 1
  ATSS.Stage1.case1 <- c(1.000, 3.000, 20.000, 35.000, 0.080,
                         0.901, 23.962, 0.736)
  expect_equal(as.numeric(ATSS_Design_Stage1(0.05, 0.20, 20, 0.10, 0.10)),
               ATSS.Stage1.case1)

  # Case 1 for Adaptive Threshold and Sample Size Simon Design (ATSS Simon) in Stage 1
  ATSS.Stage1.case2 <- c(2.000, 4.000, 18.000, 26.000,
                         0.099, 0.904, 20.130, 0.734)
  expect_equal(as.numeric(ATSS_Design_Stage1(0.10, 0.30, 18, 0.10, 0.10)),
               ATSS.Stage1.case2)
})
