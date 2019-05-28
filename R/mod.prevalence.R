
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  diag.stage <- dat$attr$diag.stage
  race <- dat$attr$race
  age <- dat$attr$age

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # Pop Size / Demog
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == 1, na.rm = TRUE)
  dat$epi$num.H[at] <- sum(race == 2, na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == 3, na.rm = TRUE)
  dat$epi$age.mean[at] <- mean(age, na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)

  # Prev / Incid
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- sum(race == 1 & status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.H[at] <- sum(race == 2 & status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.W[at] <- sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$i.prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$i.prev.dx.B[at] <- sum(race == 1 & diag.status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.dx.H[at] <- sum(race == 2 & diag.status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.dx.W[at] <- sum(race == 3 & diag.status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200

  # dat$epi$R0.mean.cs[at] <- mean(dat$attr$count.trans[status == 1], na.rm = TRUE)
  # dat$epi$R0.mean.cens[at] <- suppressWarnings(mean(tail(dat$temp$R0, 500), na.rm = TRUE))

  # Care continuum stats (primary)
  dat$epi$cc.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$cc.dx.B[at] <- sum(diag.status == 1 & race == 1, na.rm = TRUE) /
                         sum(status == 1 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.H[at] <- sum(diag.status == 1 & race == 2, na.rm = TRUE) /
                         sum(status == 1 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.W[at] <- sum(diag.status == 1 & race == 3, na.rm = TRUE) /
                         sum(status == 1 & race == 3, na.rm = TRUE)

  dat$epi$cc.dx.acute[at] <- sum(diag.status == 1 & diag.stage %in% 1:2, na.rm = TRUE) /
                             sum(diag.status == 1, na.rm = TRUE)
  dat$epi$cc.dx.chronic[at] <- sum(diag.status == 1 & diag.stage == 3, na.rm = TRUE) /
                               sum(diag.status == 1, na.rm = TRUE)

  dat$epi$cc.dx.aids[at] <- sum(diag.status == 1 & diag.stage == 4, na.rm = TRUE) /
                            sum(diag.status == 1, na.rm = TRUE)
  dat$epi$cc.dx.aids.B[at] <- sum(diag.status == 1 & diag.stage == 4 & race == 1, na.rm = TRUE) /
                              sum(diag.status == 1 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.aids.H[at] <- sum(diag.status == 1 & diag.stage == 4 & race == 2, na.rm = TRUE) /
                              sum(diag.status == 1 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.W[at] <- sum(diag.status == 1 & diag.stage == 4 & race == 3, na.rm = TRUE) /
                              sum(diag.status == 1 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4, na.rm = TRUE) /
                             sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.B[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4 & race == 1, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.H[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4 & race == 2, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  dat$epi$cc.linked1m.W[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4 & race == 3, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  dat$epi$cc.tx.any1y[at] <- sum((at - dat$attr$tx.period.last <= 52), na.rm = TRUE) /
                             sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.tx.any1y.B[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 1, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  dat$epi$cc.tx.any1y.H[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 2, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  dat$epi$cc.tx.any1y.W[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 3, na.rm = TRUE) /
                               sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$diag.status == 1, na.rm = TRUE) /
                          sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.B[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$diag.status == 1 &
                                  race == 1, na.rm = TRUE) /
                            sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.H[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$diag.status == 1 &
                                  race == 2, na.rm = TRUE) /
                            sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.W[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$diag.status == 1 &
                                  race == 3, na.rm = TRUE) /
                            sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.dur1y[at] <- 1 - (sum((at - dat$attr$vl.last.usupp) <= 52 &
                                         dat$attr$diag.status == 1, na.rm = TRUE) /
                                     sum(dat$attr$diag.status == 1, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.B[at] <- 1 - (sum((at - dat$attr$vl.last.usupp) <= 52 &
                                         dat$attr$diag.status == 1 & race == 1, na.rm = TRUE) /
                                     sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.H[at] <- 1 - (sum((at - dat$attr$vl.last.usupp) <= 52 &
                                           dat$attr$diag.status == 1 & race == 2, na.rm = TRUE) /
                                       sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.W[at] <- 1 - (sum((at - dat$attr$vl.last.usupp) <= 52 &
                                           dat$attr$diag.status == 1 & race == 3, na.rm = TRUE) /
                                       sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE))

  dat$epi$cc.HIV.mr[at] <- dat$epi$dep.HIV[at]/dat$epi$i.num[at]

  # Care continuum stats (secondary)
  dat$epi$cc.test.int[at] <- mean(at - dat$attr$last.neg.test & diag.status == 0, na.rm = TRUE)

  # dat$epi$cc.dx.delay[at] <- mean(dat$attr$diag.time - dat$attr$inf.time, na.rm = TRUE)
  # dat$epi$cc.testpy[at] <- 1-sum((at - dat$attr$last.neg.test) > 52 & status == 0,
  #     is.na(dat$attr$last.neg.test) & status == 0, na.rm = TRUE) /
  #   sum(status == 0)
  # dat$epi$cc.linked[at] <- sum(dat$attr$cuml.time.on.tx > 0, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx[at] <- sum(dat$attr$tx.status == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.ret3m[at] <- sum((at - dat$attr$tx.period.last) <= 52 &
  #       (dat$attr$tx.period.last - dat$attr$tx.period.first) > 13, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt1[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt2[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt3[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 3, na.rm = TRUE)


  # HIV stage
  dat$epi$hstage.acute[at] <- sum(dat$attr$stage %in% 1:2, na.rm = TRUE) /
                              sum(status == 1, na.rm = TRUE)
  dat$epi$hstage.chronic[at] <- sum(dat$attr$stage == 3, na.rm = TRUE) /
                                sum(status == 1, na.rm = TRUE)
  dat$epi$hstage.aids[at] <- sum(dat$attr$stage == 4, na.rm = TRUE) /
                             sum(status == 1, na.rm = TRUE)

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepCurr.hadr[at] <- sum(prepStat == 1 & prepClass == 3, na.rm = TRUE)

  # STIs
  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  ir100.rgc <- (dat$epi$incid.rgc[at]/sum(rGC == 0, dat$epi$incid.rgc[at], na.rm = TRUE))*5200
  ir100.ugc <- (dat$epi$incid.ugc[at]/sum(uGC == 0, dat$epi$incid.ugc[at], na.rm = TRUE))*5200
  dat$epi$ir100.gc[at] <- ir100.rgc + ir100.ugc
  ir100.rct <- (dat$epi$incid.rct[at]/sum(rCT == 0, dat$epi$incid.rct[at], na.rm = TRUE))*5200
  ir100.uct <- (dat$epi$incid.uct[at]/sum(uCT == 0, dat$epi$incid.uct[at], na.rm = TRUE))*5200
  dat$epi$ir100.ct[at] <- ir100.rct + ir100.uct
  dat$epi$ir100.sti[at] <- dat$epi$ir100.gc[at] + dat$epi$ir100.ct[at]

  return(dat)
}

#' @export
#' @rdname prevalence_msm
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
