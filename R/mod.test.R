
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports memoryless HIV testing for stochastic and
#' geometrically-distributed waiting times to test (constant hazard).
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
test_msm <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  hiv.test.byG.int <- dat$param$hiv.test.byG.int
  twind.int <- dat$param$test.window.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  ## Process
  elig.B <- which(race == "B" &
                  tt.traj != 1 &
                  (diag.status == 0 | is.na(diag.status)) &
                  prepStat == 0)
  rates.B <- rep(1/hiv.test.byG.int[1], length(elig.B))
  tst.B <- elig.B[rbinom(length(elig.B), 1, rates.B) == 1]

  elig.W <- which(race == "W" &
                  tt.traj != 1 &
                  (diag.status == 0 | is.na(diag.status)) &
                  prepStat == 0)
  rates.W <- rep(1/hiv.test.byG.int[2], length(elig.W))
  tst.W <- elig.W[rbinom(length(elig.W), 1, rates.W) == 1]
  tst.nprep <- c(tst.B, tst.W)


  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                    prepStat == 1 &
                    tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.prep)

  tst.pos <- tst.all[status[tst.all] == 1 & inf.time[tst.all] <= at - twind.int]
  tst.neg <- setdiff(tst.all, tst.pos)

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  return(dat)
}


#' @export
#' @rdname test_msm
dx_het <- function(dat, at) {

  # Variables
  status <- dat$attr$status
  txCD4min <- dat$attr$txCD4min
  cd4Count <- dat$attr$cd4Count
  dxStat <- dat$attr$dxStat

  # Process
  tested <- which(status == 1 & dxStat == 0 & cd4Count <= txCD4min)


  # Results
  if (length(tested) > 0) {
    dat$attr$dxStat[tested] <- 1
    dat$attr$txStat[tested] <- 0
    dat$attr$dxTime[tested] <- at
  }

  return(dat)
}

