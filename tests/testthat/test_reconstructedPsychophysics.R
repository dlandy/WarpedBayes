context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  stims <- -10:10
  
  expect_warning(bayesianGonzalezWu(-5:1/2, kappa=0.8,tauStimuli=400,tauCategory=100, responses=c(), mode= "simulation")
                 , "LeftBoundary .* larger than smallest stimulus")
  expect_warning(bayesianGonzalezWu(1:5/2, kappa=0.8,tauStimuli=400,tauCategory=100, responses=c(), mode= "simulation")
                 , "RightBoundary .* smaller than largest stimulus")
  
})

