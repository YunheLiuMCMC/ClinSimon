test_that("Test ATSS_Design_Stage2() get correct design output", {
  # Case 1 for Adaptive Threshold and Sample Size Simon Design (ATSS Simon) in Stage 2
  ATSS.Stage2.case1 <- c(1.000, 3.000, 20.000, 33.000, 0.070,
                         0.888, 23.434, 0.736)
  expect_equal(as.numeric(ATSS_Design_Stage2(0.05, 0.20, 1, 20, 33, 0.10)),
               ATSS.Stage2.case1)

  # Case 2 for Adaptive Threshold and Sample Size Simon Design (ATSS Simon) in Stage 2
  ATSS.Stage2.case2 <- c(2.000, 4.000, 18.000, 24.000,
                         0.080, 0.876, 19.597, 0.734)
  expect_equal(as.numeric(ATSS_Design_Stage2(0.10, 0.30, 2, 18, 24, 0.10)),
               ATSS.Stage2.case2)
})
