
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
  ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx
  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  # Attributes
  race <- dat$attr$race

  ## Symptomatic GC Treatment ##
  idsRGC_tx_sympt <- which(dat$attr$rGC == 1 &
                             dat$attr$rGC.infTime < at &
                             dat$attr$rGC.sympt == 1 &
                             is.na(dat$attr$rGC.tx))
  idsUGC_tx_sympt <- which(dat$attr$uGC == 1 &
                             dat$attr$uGC.infTime < at &
                             dat$attr$uGC.sympt == 1 &
                             is.na(dat$attr$uGC.tx))

  # Subset by race
  idsRGC_tx_sympt_B <- intersect(idsRGC_tx_sympt, which(race == 0))
  idsRGC_tx_sympt_W <- intersect(idsRGC_tx_sympt, which(race == 1))
  idsUGC_tx_sympt_B <- intersect(idsUGC_tx_sympt, which(race == 0))
  idsUGC_tx_sympt_W <- intersect(idsUGC_tx_sympt, which(race == 1))

  # Collect over site
  idsGC_tx_sympt_B <- union(idsRGC_tx_sympt_B, idsUGC_tx_sympt_B)
  idsGC_tx_sympt_W <- union(idsRGC_tx_sympt_W, idsUGC_tx_sympt_W)

  # Treatment by race
  txGC_sympt_B <- idsGC_tx_sympt_B[which(rbinom(length(idsGC_tx_sympt_B), 1,
                                                gc.sympt.prob.tx[1]) == 1)]
  txGC_sympt_W <- idsGC_tx_sympt_W[which(rbinom(length(idsGC_tx_sympt_W), 1,
                                                gc.sympt.prob.tx[2]) == 1)]
  txGC_sympt <- union(txGC_sympt_B, txGC_sympt_W)

  # Subset by site
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)

  ## Asymptomatic GC Treatment ##
  idsRGC_tx_asympt <- which(dat$attr$rGC == 1 &
                              dat$attr$rGC.infTime < at &
                              dat$attr$rGC.sympt == 0 &
                              is.na(dat$attr$rGC.tx) &
                              dat$attr$prepStat == 0)
  idsUGC_tx_asympt <- which(dat$attr$uGC == 1 &
                              dat$attr$uGC.infTime < at &
                              dat$attr$uGC.sympt == 0 &
                              is.na(dat$attr$uGC.tx) &
                              dat$attr$prepStat == 0)

  # Subset by race
  idsRGC_tx_asympt_B <- intersect(idsRGC_tx_asympt, which(race == 0))
  idsRGC_tx_asympt_W <- intersect(idsRGC_tx_asympt, which(race == 1))
  idsUGC_tx_asympt_B <- intersect(idsUGC_tx_asympt, which(race == 0))
  idsUGC_tx_asympt_W <- intersect(idsUGC_tx_asympt, which(race == 1))

  # Collect over site
  idsGC_tx_asympt_B <- union(idsRGC_tx_asympt_B, idsUGC_tx_asympt_B)
  idsGC_tx_asympt_W <- union(idsRGC_tx_asympt_W, idsUGC_tx_asympt_W)

  # Treatment by race
  txGC_asympt_B <- idsGC_tx_asympt_B[which(rbinom(length(idsGC_tx_asympt_B), 1,
                                                  gc.asympt.prob.tx[1]) == 1)]
  txGC_asympt_W <- idsGC_tx_asympt_W[which(rbinom(length(idsGC_tx_asympt_W), 1,
                                                  gc.asympt.prob.tx[2]) == 1)]
  txGC_asympt <- union(txGC_asympt_B, txGC_asympt_W)

  # Subset by site
  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)

  ## All Treated GC ##

  # IDs of men sucessfully treated
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  # IDs of men eligible for treatment
  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  ## Symptomatic CT Treatment ##
  idsRCT_tx_sympt <- which(dat$attr$rCT == 1 &
                             dat$attr$rCT.infTime < at &
                             dat$attr$rCT.sympt == 1 &
                             is.na(dat$attr$rCT.tx))
  idsUCT_tx_sympt <- which(dat$attr$uCT == 1 &
                             dat$attr$uCT.infTime < at &
                             dat$attr$uCT.sympt == 1 &
                             is.na(dat$attr$uCT.tx))

  # Subset by race
  idsRCT_tx_sympt_B <- intersect(idsRCT_tx_sympt, which(race == 0))
  idsRCT_tx_sympt_W <- intersect(idsRCT_tx_sympt, which(race == 1))
  idsUCT_tx_sympt_B <- intersect(idsUCT_tx_sympt, which(race == 0))
  idsUCT_tx_sympt_W <- intersect(idsUCT_tx_sympt, which(race == 1))

  # Collect over site
  idsCT_tx_sympt_B <- union(idsRCT_tx_sympt_B, idsUCT_tx_sympt_B)
  idsCT_tx_sympt_W <- union(idsRCT_tx_sympt_W, idsUCT_tx_sympt_W)

  # Treatment by race
  txCT_sympt_B <- idsCT_tx_sympt_B[which(rbinom(length(idsCT_tx_sympt_B), 1,
                                                ct.sympt.prob.tx[1]) == 1)]
  txCT_sympt_W <- idsCT_tx_sympt_W[which(rbinom(length(idsCT_tx_sympt_W), 1,
                                                ct.sympt.prob.tx[2]) == 1)]
  txCT_sympt <- union(txCT_sympt_B, txCT_sympt_W)

  # Subset by site
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)


  ## Asymptomatic CT Treatment ##
  idsRCT_tx_asympt <- which(dat$attr$rCT == 1 &
                              dat$attr$rCT.infTime < at &
                              dat$attr$rCT.sympt == 0 &
                              is.na(dat$attr$rCT.tx) &
                              dat$attr$prepStat == 0)
  idsUCT_tx_asympt <- which(dat$attr$uCT == 1 &
                              dat$attr$uCT.infTime < at &
                              dat$attr$uCT.sympt == 0 &
                              is.na(dat$attr$uCT.tx) &
                              dat$attr$prepStat == 0)

  # Subset by race
  idsRCT_tx_asympt_B <- intersect(idsRCT_tx_asympt, which(race == 0))
  idsRCT_tx_asympt_W <- intersect(idsRCT_tx_asympt, which(race == 1))
  idsUCT_tx_asympt_B <- intersect(idsUCT_tx_asympt, which(race == 0))
  idsUCT_tx_asympt_W <- intersect(idsUCT_tx_asympt, which(race == 1))

  # Collect over site
  idsCT_tx_asympt_B <- union(idsRCT_tx_asympt_B, idsUCT_tx_asympt_B)
  idsCT_tx_asympt_W <- union(idsRCT_tx_asympt_W, idsUCT_tx_asympt_W)

  # Treatment by race
  txCT_asympt_B <- idsCT_tx_asympt_B[which(rbinom(length(idsCT_tx_asympt_B), 1,
                                                  ct.asympt.prob.tx[1]) == 1)]
  txCT_asympt_W <- idsCT_tx_asympt_W[which(rbinom(length(idsCT_tx_asympt_W), 1,
                                                  ct.asympt.prob.tx[2]) == 1)]
  txCT_asympt <- union(txCT_asympt_B, txCT_asympt_W)

  # Subset by site
  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)

  ## All Treated CT ##
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)


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
  idsRCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rCT == 1 &
                                      dat$attr$rCT.infTime < at &
                                      is.na(dat$attr$rCT.tx.prep)))
  idsUCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uCT == 1 &
                                      dat$attr$uCT.infTime < at &
                                      is.na(dat$attr$uCT.tx.prep)))

  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]


  ## Update Attributes ##
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1

  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1

  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1

  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1

  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1

  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1

  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1


  ## Add tx at other anatomical site ##
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) &
                          dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) &
                          dat$attr$uGC == 1)] <- 1

  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1) &
                          dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1) &
                          dat$attr$uCT == 1)] <- 1

  return(dat)
}
