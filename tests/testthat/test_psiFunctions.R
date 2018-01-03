context("Psychophysical Transformation Functions")

test_that("psiLinear Shifts and scales", {
  stims <- -10:10
  
  expect_identical(psiLinear(stims, shift=10, scaling=2), (stims)*2+10)
  expect_equal(psiLinear(stims), stims)
})

test_that("psiIdentical is", {
  stims <- -10:10
  
  expect_equal(psiIdentity(stims), stims)
})