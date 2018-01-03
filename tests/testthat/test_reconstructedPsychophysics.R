context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  stims <- -10:10
  
  expect_warning(bayesianGonzalezWu(-5:1, kappa=0.8,tauStimuli=400,tauCategory=100, responses=c(), mode= "simulation")
                 , "LeftBoundary (-1e-10) larger than smallest stimulus (-2.5)")
  
  
})

