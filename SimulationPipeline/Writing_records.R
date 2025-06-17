cat(" Writing records \n")
# NOTE: Training records are collected with a sliding-window process.
# Cumulation of the records starts in the burn-in year 'startTrainPop'.
# Once the burn-in period is over, the sliding-window process removes the oldest
# records.


if(year >= startTrainPop){

  ######------------------------------------------------
  ##>>>>>  Pull population
  ######------------------------------------------------

  if (year == 12) {
    DHTrain1@fixEff <- as.integer(rep(paste0(year,1L), nInd(DHTrain1)))
    DHTrain2@fixEff <- as.integer(rep(paste0(year,2L), nInd(DHTrain2)))
    DHTrain3@fixEff <- as.integer(rep(paste0(year,3L), nInd(DHTrain3)))
    
    trainPop = c(DHTrain1, DHTrain2, DHTrain3)
    

  } else if (year == 13 | year == 14 | year == 15) {

    DHTrain1@fixEff <- as.integer(rep(paste0(year,1L), nInd(DHTrain1)))
    DHTrain2@fixEff <- as.integer(rep(paste0(year,2L), nInd(DHTrain2)))
    DHTrain3@fixEff <- as.integer(rep(paste0(year,3L), nInd(DHTrain3)))
    trainYear = c(DHTrain1, DHTrain2, DHTrain3)
    
    trainPop = c(trainPop, trainYear)

  } else {

    DHTrain1@fixEff <- as.integer(rep(paste0(year,1L), nInd(DHTrain1)))
    DHTrain2@fixEff <- as.integer(rep(paste0(year,2L), nInd(DHTrain2)))
    DHTrain3@fixEff <- as.integer(rep(paste0(year,3L), nInd(DHTrain3)))
    trainYear = c(DHTrain1, DHTrain2, DHTrain3)
    

    trainPop <- trainPop[-c(1:nInd(trainYear))]
    trainPop = c(trainPop, trainYear)

  }

}

