context("bayesianReconstructionFunctions")

test_that("bayesianGonzalezWu gives warnings and return values", {
  a <- with(boundedProportionSimulatedData,
            fitWarpedBayesModel(bayesianGonzalezWu, 
                              stimulus
                           ,  response
                           , initialPars = c(kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
                           , responseGrid = sort(unique(response))
                           , fixedPars=c() ))
  expect_lt(a$kappa[1], 0.3)
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





test_that("Guess Model returns reasonable values", {
  expect_lt(max(abs((guessingModel(bayesianGonzalezWu))(1:9/10, mode="prediction", guessRate=0.5,  responseGrid=1:9/10)
                    - c(0.3000045, 0.4097168, 0.4464102, 0.4745967, 0.5000000, 0.5254033, 0.5535898, 0.5902832, 0.6999955))),
            1e-4)
  expect_lt(max(abs((guessingModel(bayesianGonzalezWu))(1:9/10, mode="prediction", guessRate=0.1,  responseGrid=1:9/10)
                    - c(0.1400080, 0.3374902, 0.4035383, 0.4542740, 0.5000000, 0.5457260, 0.5964617, 0.6625098, 0.8599920))),
            1e-4)
  expect_lt(max(abs((guessingModel(bayesianGonzalezWu))(1:9/10, mode="prediction", guessRate=0.0,  responseGrid=1:9/10)
                    - c(0.1000089, 0.3194335, 0.3928203, 0.4491933, 0.5000000, 0.5508067, 0.6071797, 0.6805665, 0.8999911))),
            1e-4)
  expect_lt(max(abs((guessingModel(bayesianGonzalezWu))(1:9/10, mode="prediction", guessRate=1.0,  responseGrid=1:9/10)
                    - c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))),
            1e-4)
  
  
  a <- with(boundedProportionWithGuessingSimulatedData,
            fitWarpedBayesModel(guessingModel(bayesianGonzalezWu), 
                                stimulus
                                ,  response
                                , initialPars = c(guessRate = 0, kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
                                , responseGrid = sort(unique(response))
                                , fixedPars=c() ))
  
  expect_lt(a$kappa[1], 0.30)
  expect_gt(a$kappa[1], 0.10)
  expect_lt(a$tauStimuli[1], 450)
  expect_gt(a$tauStimuli[1], 350)
  expect_lt(a$tauCategory[1], 110)
  expect_gt(a$tauCategory[1], 80)
  expect_lt(a$guessRate[1], .40)
  expect_gt(a$guessRate[1], .25)
})






