context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  a <- with(boundedProportionSimulatedData,
            fitWarpedBayesModel(bayesianGonzalezWu, 
                              stimulus
                           ,  response
                           , initialPars = c(kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
                           , responseGrid = sort(unique(response))
                           , fixedPars=c() ))
  expect_lt(a$kappa[1], 0.2)
  expect_gt(a$kappa[1], 0.1)
  expect_lt(a$tauStimuli[1], 450)
  expect_gt(a$tauStimuli[1], 350)
  expect_lt(a$tauCategory[1], 120)
  expect_gt(a$tauCategory[1], 80)
  
})




test_that("bayesianGonzalezWu gives reasonable values", {
  fakeStims <- 1:10/10
  fakeData <- bayesianGonzalezWu(fakeStims, kappa= 0,tauStimuli=786,tauCategory=262, rightBoundaryExpansion=-8, mode="prediction")
  expect_lt(max(abs(fakeData- c(0.1000000, 0.2563537, 0.3529170, 0.4356332, 0.5124735, 0.5876090, 0.6644375, 0.7471244, 0.8436163, 0.9979267))),
            1e-4)
})



test_that("bayesianSpatialMemoryLandyCrawfordCorbin2017 gives warnings and return values", {
  a <- with(spatialMemorySimulatedData,
            fitWarpedBayesModel(bayesianSpatialMemoryLandyCrawfordCorbin2017, 
                                stimulus
                                ,  response
                                , initialPars = c(kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
                                , responseGrid = sort(unique(response))
                                , fixedPars=c() ))
  expect_lt(a$kappa[1], 0.55)
  expect_gt(a$kappa[1], 0.45)
  expect_lt(a$tauStimuli[1], 90)
  expect_gt(a$tauStimuli[1], 65)
  expect_lt(a$tauCategory[1], 30)
  expect_gt(a$tauCategory[1], 15)
  
})
