context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  stims <- -10:10
  
  #expect_warning(bayesianGonzalezWu(-5:1/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode = "simulation")
   #              , "leftBoundaryObjective .* larger than smallest stimulus")
  #expect_warning(bayesianGonzalezWu(1:5/2, kappa=0.8,tauStimuli=400,tauCategory=100, mode= "simulation")
    #             , "rightBoundaryObjective .* smaller than largest stimulus")
  
})



test_that("bayesianGonzalezWu gives reasonable values", {
  fakeStims <- 1:10/10
  fakeData <- bayesianGonzalezWu(fakeStims, kappa= 0,tauStimuli=786,tauCategory=262, rightBoundaryExpansion=-8, mode="prediction")
  expect_lt(max(abs(fakeData- c(0.1000000, 0.2563537, 0.3529170, 0.4356332, 0.5124735, 0.5876090, 0.6644375, 0.7471244, 0.8436163, 0.9979267))),
            1e-4)
})

