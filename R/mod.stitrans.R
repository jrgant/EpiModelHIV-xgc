
#' @title STI Transmission Module
#'
#' @description Stochastically simulates GC transmission given the current
#'              state of the edgelist.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stitrans_msm <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  # Acquisition probabilities given contact with infected man
  rgc.tprob <- dat$param$rgc.tprob
  ugc.tprob <- dat$param$ugc.tprob
  pgc.tprob <- dat$param$pgc.tprob

  # Probability of symptoms given infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  pgc.sympt.prob <- dat$param$pgc.sympt.prob

  # Relative risk of infection given condom use during act
  sti.cond.eff <- dat$param$sti.cond.eff
  sti.cond.fail <- dat$param$sti.cond.fail


  # Attributes ----------------------------------------------------------

  race <- dat$attr$race

  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  pGC <- dat$attr$pGC

  # n Times infected
  rGC.timesInf <- dat$attr$rGC.timesInf
  uGC.timesInf <- dat$attr$uGC.timesInf
  pGC.timesInf <- dat$attr$pGC.timeInf

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  pGC.infTime <- dat$attr$pGC.infTime

  # Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  pGC.sympt <- dat$attr$pGC.sympt


  # Pull act list
  al <- dat$temp$al

  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive


  # Rectal GC -----------------------------------------------------------

  # Requires: uGC in insertive man, and no rGC in receptive man
  p1Inf_rgc <- which(uGC[al[, "p1"]] == 1 & uGC.infTime[al[, "p1"]] < at &
                     rGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rgc <- which(uGC[al[, "p2"]] == 1 & uGC.infTime[al[, "p2"]] < at &
                     rGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rgc <- c(p1Inf_rgc, p2Inf_rgc)

  # UAI modifier
  uai_rgc <- al[allActs_rgc, "uai"]
  tprob_rgc <- rep(rgc.tprob, length(allActs_rgc))

  # Transform to log odds
  tlo_rgc <- log(tprob_rgc/(1 - tprob_rgc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_rgc, "p1"]], race[al[p2Inf_rgc, "p2"]])
  condom.rr <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_rgc[uai_rgc == 0] <- tlo_rgc[uai_rgc == 0] + log(condom.rr[uai_rgc == 0])

  # Back-transform to probability
  tprob_rgc <- plogis(tlo_rgc)

  # Stochastic transmission
  trans_rgc <- rbinom(length(allActs_rgc), 1, tprob_rgc)

  # Determine the infected partner
  idsInf_rgc <- NULL
  if (sum(trans_rgc) > 0) {
    transAL_rgc <- al[allActs_rgc[trans_rgc == 1], , drop = FALSE]
    idsInf_rgc <- c(intersect(al[p1Inf_rgc, "p2"], transAL_rgc[, "p2"]),
                    intersect(al[p2Inf_rgc, "p1"], transAL_rgc[, "p1"]))
    stopifnot(all(rGC[idsInf_rgc] == 0))
  }

  # Update attributes
  rGC[idsInf_rgc] <- 1
  rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- rbinom(length(idsInf_rgc), 1, rgc.sympt.prob)
  rGC.timesInf[idsInf_rgc] <- rGC.timesInf[idsInf_rgc] + 1


  # Urethral GC ---------------------------------------------------------

  # Requires: rGC in receptive man, and no uGC in insertive man
  p1Inf_ugc <- which(rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
                     uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_ugc <- which(rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
                     uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_ugc <- c(p1Inf_ugc, p2Inf_ugc)

  # UAI modifier
  uai_ugc <- al[allActs_ugc, "uai"]
  tprob_ugc <- rep(ugc.tprob, length(allActs_ugc))

  # Transform to log odds
  tlo_ugc <- log(tprob_ugc/(1 - tprob_ugc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_ugc, "p2"]], race[al[p2Inf_ugc, "p1"]])
  condom.rr <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_ugc[uai_ugc == 0] <- tlo_ugc[uai_ugc == 0] + log(condom.rr[uai_ugc == 0])

  # Back-transform to probability
  tprob_ugc <- plogis(tlo_ugc)

  # Stochastic transmission
  trans_ugc <- rbinom(length(allActs_ugc), 1, tprob_ugc)

  # Determine the newly infected partner
  idsInf_ugc <- NULL
  if (sum(trans_ugc) > 0) {
    transAL_ugc <- al[allActs_ugc[trans_ugc == 1],  , drop = FALSE]
    idsInf_ugc <- c(intersect(al[p1Inf_ugc, "p2"], transAL_ugc[, "p2"]),
                    intersect(al[p2Inf_ugc, "p1"], transAL_ugc[, "p1"]))
    stopifnot(all(uGC[idsInf_ugc] == 0))
  }

  # Update attributes
  uGC[idsInf_ugc] <- 1
  uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)
  uGC.timesInf[idsInf_ugc] <- uGC.timesInf[idsInf_ugc] + 1

  # Pharyngeal GC ----------------------------------------------------------

  # TODO 2020-04-20: Fill out pharyngeal GC section


  # Output --------------------------------------------------------------

  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$pGC <- pGC

  dat$attr$rGC.infTime <- rGC.infTime
  dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$pGC.infTime <- pGC.infTime

  dat$attr$rGC.timesInf <- rGC.timesInf
  dat$attr$uGC.timesInf <- uGC.timesInf
  dat$attr$pGC.timesInf <- pGC.timesInf

  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$pGC.sympt <- pGC.sympt

  # Summary stats

  ## GC incidence by anatomic site
  dat$epi$incid.gc[at] <- length(idsInf_rgc) + length(idsInf_ugc) + length(idsInf_pgc)
  dat$epi$incid.rgc[at] <- length(idsInf_rgc)
  dat$epi$incid.ugc[at] <- length(idsInf_ugc)
  dat$epi$incid.pgc[at] <- length(idsInf_pgc)

  ## by race/ethnicity
  dat$epi$incid.gc.B[at] <- length(
    intersect(union(idsInf_rgc, idsInf_ugc, idsInf_pgc), which(race == 1))
  )
  dat$epi$incid.gc.H[at] <- length(
    intersect(union(idsInf_rgc, idsInf_ugc, idsInf_pgc), which(race == 2))
  )
  dat$epi$incid.gc.O[at] <- length(
    intersect(union(idsInf_rgc, idsInf_ugc, idsInf_pgc), which(race == 3))
  )
  dat$epi$incid.gc.W[at] <- length(
    intersect(union(idsInf_rgc, idsInf_ugc, idsInf_pgc), which(race == 4))
  )

  # Check all infected have all STI attributes
  stopifnot(
    all(!is.na(rGC.infTime[rGC == 1])),
    all(!is.na(rGC.sympt[rGC == 1])),
    all(!is.na(uGC.infTime[uGC == 1])),
    all(!is.na(uGC.sympt[uGC == 1])),
    all(!is.na(pGC.infTime[pGC == 1])),
    all(!is.na(pGC.sympt[pGC == 1]))
  )

  return(dat)
}
