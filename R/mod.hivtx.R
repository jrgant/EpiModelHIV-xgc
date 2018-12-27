
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Persons enter into the simulation with one of four ART "patterns": never
#' tested, tested but never treated, treated and achieving partial HIV viral
#' suppression, or treated with full viral suppression (these types are stored
#' as individual-level attributes in \code{tt.traj}). This module initiates ART
#' for treatment naive persons in the latter two types, and then cycles them on
#' and off treatment conditional on empirical race-specific adherence rates. ART
#' initiation, non-adherence, and restarting are all stochastically simulated
#' based on binomial statistical models.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{tx.status},
#' \code{tx.init.time}, \code{cum.time.on.tx}, \code{cum.time.off.tx} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
hivtx_msm <- function(dat, at) {

  # Attributes
  race <- dat$attr$race
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  cum.time.off.tx <- dat$attr$cum.time.off.tx
  stage <- dat$attr$stage

  # Parameters
  tx.init.prob <- dat$param$tx.init.prob
  tx.halt.prob <- dat$param$tx.halt.prob
  tx.reinit.prob <- dat$param$tx.reinit.prob

  ## Initiation
  tx.init.elig <- which(status == 1 &
                        tx.status == 0 &
                        diag.status == 1 &
                        cum.time.on.tx == 0 &
                        stage != 4)
  rates <- tx.init.prob[race[tx.init.elig] + 1]
  tx.init <- tx.init.elig[rbinom(length(tx.init.elig), 1, rates) == 1]

  ## Halting
  tx.halt.elig <- which(tx.status == 1)
  rates <- tx.halt.prob[race[tx.halt.elig] + 1]
  tx.halt <- tx.halt.elig[rbinom(length(tx.halt.elig), 1, rates) == 1]

  ## Restarting
  tx.reinit.elig <- which(tx.status == 0 &
                          cum.time.on.tx > 0 &
                          stage != 4)
  rates <- tx.reinit.prob[race[tx.reinit.elig] + 1]
  tx.reinit <- tx.reinit.elig[rbinom(length(tx.reinit.elig), 1, rates) == 1]

  ## Update treatment time
  cum.time.on.tx[which(tx.status == 1)] <- cum.time.on.tx[which(tx.status == 1)] + 1
  cum.time.off.tx[which(tx.status == 0)] <- cum.time.off.tx[which(tx.status == 0)] + 1

  ## Update Attributes
  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.status[tx.halt] <- 0
  dat$attr$tx.status[tx.reinit] <- 1

  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx

  return(dat)
}


#' @export
#' @rdname hivtx_msm
tx_het <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  dxStat <- dat$attr$dxStat
  txStat <- dat$attr$txStat
  txStartTime <- dat$attr$txStartTime
  txStops <- dat$attr$txStops
  txTimeOn <- dat$attr$txTimeOn
  txTimeOff <- dat$attr$txTimeOff
  txCD4start <- dat$attr$txCD4start

  cd4Count <- dat$attr$cd4Count
  tx.elig.cd4 <- dat$param$tx.elig.cd4
  tx.coverage <- dat$param$tx.coverage

  txType <- dat$attr$txType
  tx.adhere.full <- dat$param$tx.adhere.full
  tx.adhere.part <- dat$param$tx.adhere.part


  # Start tx for tx naive ---------------------------------------------------

  ## Calculate tx coverage
  allElig <- which((cd4Count < tx.elig.cd4 | !is.na(txStartTime)))
  txCov <- sum(!is.na(txStartTime[allElig]))/length(allElig)
  if (is.nan(txCov)) {
    txCov <- 0
  }

  idsElig <- which(dxStat == 1 & txStat == 0 &
                   is.na(txStartTime) & cd4Count < tx.elig.cd4)
  nElig <- length(idsElig)
  idsTx <- NULL


  ## Treatment coverage
  nStart <- max(0, min(nElig, round((tx.coverage - txCov) * length(allElig))))
  if (nStart > 0) {
    idsTx <- ssample(idsElig, nStart)
  }


  ## Treatment type assignment
  if (length(idsTx) > 0) {
    needtxType <- which(is.na(txType[idsTx]))
    if (length(needtxType) > 0) {
      txType[idsTx[needtxType]] <- rbinom(length(needtxType), 1, tx.adhere.full)
    }
    if (tx.adhere.part == 0) {
      idsTx <- intersect(idsTx, which(txType == 1))
    }
  }

  if (length(idsTx) > 0) {
    txStat[idsTx] <- 1
    txStartTime[idsTx] <- at
    txStops[idsTx] <- 0
    txTimeOn[idsTx] <- 0
    txTimeOff[idsTx] <- 0
    txCD4start[idsTx] <- cd4Count[idsTx]
  }


  # Stop tx -----------------------------------------------------------------
  idsStop <- NULL
  idsEligStop <- which(dat$attr$txStat == 1 & txType == 0)
  nEligStop <- length(idsEligStop)
  if (nEligStop > 0) {
    vecStop <- which(rbinom(nEligStop, 1, (1 - tx.adhere.part)) == 1)
    if (length(vecStop) > 0) {
      idsStop <- idsEligStop[vecStop]
      txStat[idsStop] <- 0
      txStops[idsStop] <- txStops[idsStop] + 1
    }
  }


  # Restart tx --------------------------------------------------------------
  idsRest <- NULL
  idsEligRest <- which(dat$attr$txStat == 0 & txStops > 0)
  nEligRest <- length(idsEligRest)
  if (nEligRest > 0) {
    vecRes <- which(rbinom(nEligRest, 1, tx.adhere.part) == 1)
    if (length(vecRes) > 0) {
      idsRest <- idsEligRest[vecRes]
      txStat[idsRest] <- 1
      dat$attr$vlSlope[idsRest] <- NA
    }
  }


  # Output ------------------------------------------------------------------
  idsOnTx <- which(txStat == 1)
  idsOffTx <- which(txStat == 0 & !is.na(txStartTime))
  txTimeOn[idsOnTx] <- txTimeOn[idsOnTx] + 1
  txTimeOff[idsOffTx] <- txTimeOff[idsOffTx] + 1

  dat$attr$txStat <- txStat
  dat$attr$txStartTime <- txStartTime
  dat$attr$txStops <- txStops
  dat$attr$txTimeOn <- txTimeOn
  dat$attr$txTimeOff <- txTimeOff
  dat$attr$txType <- txType
  dat$attr$txCD4start <- txCD4start

  return(dat)
}

