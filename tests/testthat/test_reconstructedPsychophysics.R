context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  stims <- -10:10
  
  expect_warning(bayesianGonzalezWu(-5:1/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode = "simulation")
                 , "LeftBoundaryObj .* larger than smallest stimulus")
  expect_warning(bayesianGonzalezWu(1:5/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode= "simulation")
                 , "RightBoundaryObj .* smaller than largest stimulus")
  
})



test_that("bayesianGonzalezWu gives reasonable values", {
  fakeStims <- 1:10/10
  fakeData <- bayesianGonzalezWu(fakeStims, kappa= -0.5,tauStimuli=786,tauCategory=262, rightBoundaryExpansion=-8, mode="prediction")
  expect_lt(max(abs(fakeData- c(0.1001413, 0.2408544, 0.3308149, 0.4097610, 0.4847134, 0.5595559, 0.6377240, 0.7237672, 0.8267800, 0.997780))),
            1e-4)
})

