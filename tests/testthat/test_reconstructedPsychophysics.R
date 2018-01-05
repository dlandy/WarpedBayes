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
  expect_lt(max(abs(fakeData- c(0.1631952, 0.2610112, 0.3453110, 0.4232860, 0.4986531, 0.5741620, 0.6526095, 0.7379117, 0.8380115, 0.9978551))),
            1e-4)
})

