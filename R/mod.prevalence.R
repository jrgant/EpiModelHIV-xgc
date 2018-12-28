
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
  race <- dat$attr$race

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == 0, na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == 1, na.rm = TRUE)

  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)

  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.prep[at] <- sum(status == 1 & (prepStat == 1), na.rm = TRUE) /
    sum(prepStat == 1, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, na.rm = TRUE)) * 5200

  # Care continuum stats
  dat$epi$cc.dx[at] <- sum(dat$attr$diag.status == 1, na.rm = TRUE) /
    sum(status == 1, na.rm = TRUE)
  dat$epi$cc.linked[at] <- sum(dat$attr$cum.time.on.tx > 0, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.linked1m[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.tx[at] <- sum(dat$attr$tx.status == 1, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.tx.any1y[at] <- sum((at - dat$attr$tx.period.last <= 52), na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.tx.ret3m[at] <- sum((at - dat$attr$tx.period.last) <= 52 &
        (dat$attr$tx.period.last - dat$attr$tx.period.first) > 13, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.vsupp[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$diag.status == 1, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.tt1[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 1, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.tt2[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 2, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.tt3[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 3, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 3, na.rm = TRUE)
  dat$epi$cc.vsupp.dur1y[at] 1-(sum((at - dat$attr$vl.last.usupp) <= 52 &
                                      dat$attr$diag.status == 1, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1, na.rm = TRUE))


  # HIV stage
  dat$epi$hstage.acute[at] <- sum(dat$attr$stage %in% 1:2, na.rm = TRUE) /
    sum(status == 1, na.rm = TRUE)
  dat$epi$hstage.chronic[at] <- sum(dat$attr$stage == 3, na.rm = TRUE) /
    sum(status == 1, na.rm = TRUE)
  dat$epi$hstage.aids[at] <- sum(dat$attr$stage == 4, na.rm = TRUE) /
    sum(status == 1, na.rm = TRUE)

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)

  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$ir100.gc[at] <- (dat$epi$incid.gc[at] /
                             (sum(rGC == 0, na.rm = TRUE) +
                                sum(uGC == 0, na.rm = TRUE))) * 5200
  dat$epi$ir100.ct[at] <- (dat$epi$incid.ct[at] /
                             (sum(rCT == 0, na.rm = TRUE) +
                                sum(uCT == 0, na.rm = TRUE))) * 5200


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
