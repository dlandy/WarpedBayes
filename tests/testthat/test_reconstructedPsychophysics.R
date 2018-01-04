context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  stims <- -10:10
  
  expect_warning(bayesianGonzalezWu(-5:1/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode = "simulation")
                 , "leftBoundaryObjective .* larger than smallest stimulus")
  expect_warning(bayesianGonzalezWu(1:5/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode= "simulation")
                 , "rightBoundaryObjective .* smaller than largest stimulus")
  
})



test_that("bayesianGonzalezWu gives reasonable values", {
  fakeStims <- 1:10/10
  fakeData <- bayesianGonzalezWu(fakeStims, kappa= 0,tauStimuli=786,tauCategory=262, rightBoundaryExpansion=-8, mode="prediction")
  expect_lt(max(abs(fakeData- c(0.1614052, 0.2612281, 0.3463018, 0.4245918, 0.5000377, 0.5754796, 0.6537560, 0.7388007, 0.8385564, 0.9980160))),
            1e-4)
})

