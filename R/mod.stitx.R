
#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stitx_msm <- function(dat, at) {

  # Parameters
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ## ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx
  ## ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  # Attributes
  race <- dat$attr$race

  ## Symptomatic GC Treatment ##
  idsRGC_tx_sympt <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    dat$attr$rGC.sympt == 1 &
    is.na(dat$attr$rGC.tx)
  )

  idsUGC_tx_sympt <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    dat$attr$uGC.sympt == 1 &
    is.na(dat$attr$uGC.tx)
  )

  idsPGC_tx_sympt <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    dat$attr$pGC.sympt == 1 &
    is.na(dat$attr$uGC.tx)
  )

  # Subset by race
  idsGC_tx_sympt <- union(idsRGC_tx_sympt, idsUGC_tx_sympt, idsPGC_tx_sympt)
  races <- sort(unique(race[idsGC_tx_sympt]))
  txGC_sympt <- rep(NA, length(idsGC_tx_sympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_sympt] == i)
    txGC_sympt[ids.race] <- rbinom(length(ids.race), 1, gc.sympt.prob.tx[i])
  }
  ids_txGC_sympt <- idsGC_tx_sympt[which(txGC_sympt == 1)]

  # Subset by anatomic site
  txRGC_sympt <- intersect(idsRGC_tx_sympt, ids_txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, ids_txGC_sympt)
  txPGC_sympt <- intersect(idsPGC_tx_sympt, ids_txGC_sympt)

  ## Asymptomatic GC Treatment ##
  idsRGC_tx_asympt <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    dat$attr$rGC.sympt == 0 &
    is.na(dat$attr$rGC.tx) &
    dat$attr$prepStat == 0
  )

  idsUGC_tx_asympt <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    dat$attr$uGC.sympt == 0 &
    is.na(dat$attr$uGC.tx) &
    dat$attr$prepStat == 0
  )

  idsPGC_tx_asympt <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    dat$attr$pGC.sympt == 0 &
    is.na(dat$attr$pGC.tx) &
    dat$attr$prepStat == 0
  )

  # Subset by race
  idsGC_tx_asympt <- union(idsRGC_tx_asympt, idsUGC_tx_asympt, idsPGC_tx_asympt)
  races <- sort(unique(race[idsGC_tx_asympt]))
  txGC_asympt <- rep(NA, length(idsGC_tx_asympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_asympt] == i)
    txGC_asympt[ids.race] <- rbinom(length(ids.race), 1, gc.asympt.prob.tx[i])
  }
  ids_txGC_asympt <- idsGC_tx_asympt[which(txGC_asympt == 1)]

  # Subset by site
  txRGC_asympt <- intersect(idsRGC_tx_asympt, ids_txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, ids_txGC_asympt)
  txPGC_asympt <- intersect(idsPGC_tx_asympt, ids_txGC_asympt)

  ## All Treated GC ##

  # IDs of men sucessfully treated
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)
  txPGC <- union(txPGC_sympt, txPGC_asympt)

  # IDs of men eligible for treatment
  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)
  idsPGC_tx <- union(idsPGC_tx_sympt, idsPGC_tx_asympt)

  ## Interval-based treatment for MSM on PrEP ##
  idsSTI_screen <- which(dat$attr$prepStartTime == at |
                           (at - dat$attr$prepLastStiScreen >= prep.sti.screen.int))

  dat$attr$prepLastStiScreen[idsSTI_screen] <- at


  idsRGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rGC == 1 &
                                    dat$attr$rGC.infTime < at &
                                    is.na(dat$attr$rGC.tx.prep)))
  idsUGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uGC == 1 &
                                    dat$attr$uGC.infTime < at &
                                    is.na(dat$attr$uGC.tx.prep)))

  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]

  ## Update Attributes ##
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1

  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1

  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1


  ## Add tx at other anatomical site ##
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) &
                          dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) &
                          dat$attr$uGC == 1)] <- 1

  return(dat)
}
