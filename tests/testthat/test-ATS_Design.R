test_that("Test ATS_Design() get correct design output", {
  # Case 1 for Adaptive Threshold Simon Design (ATS Simon)
  ATS.case1 <- c(3.000, 10.000, 17.000, 34.000, 0.091,
                0.059, 0.847, 24.669, 0.549)
  expect_equal(as.numeric(ATS_Design(19, 36, 17, 34, 3, 10, 0.20, 0.40, 0.1)),
               ATS.case1)

  # Case 2 for Adaptive Threshold Simon Design (ATS Simon)
  ATS.case2 <- c(2.000, 14.000, 11.000, 41.000, 0.088,
                      0.060, 0.854, 27.344, 0.455)

  expect_equal(as.numeric(ATS_Design(14, 44, 11, 41, 3, 14, 0.25, 0.45, 0.1)),
               ATS.case2)
})
