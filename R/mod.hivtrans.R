
# MSM -----------------------------------------------------------------

#' @title HIV transmission Module
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This is the final substantive HIV-related function that occurs within
#' the time loop at each time step. This function takes the discordant
#' edgelist and calculates a transmission probability for each row
#' (one sexual act) between dyads on the network. After transmission
#' events, individual-level attributes for the infected persons are
#' updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' insertive versus receptive role, viral load of the infected partner, an
#' acute stage infection excess risk, and condom use.
#'
#' Given these transmission probabilities, transmission is stochastically
#' simulating by drawing from a binomial distribution for each act conditional
#' on the per-act probability.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
hivtrans_msm <- function(dat, at) {

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  circ <- dat$attr$circ
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  race <- dat$attr$race
  age.grp <- dat$attr$age.grp

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  trans.scale <- dat$param$trans.scale
  acute.rr <- dat$param$acute.rr

  cond.eff <- dat$param$cond.eff
  cond.fail <- dat$param$cond.fail

  circ.rr <- dat$param$circ.rr
  prep.hr <- dat$param$prep.adhr.hr
  hiv.ugc.rr <- dat$param$hiv.ugc.rr
  hiv.rgc.rr <- dat$param$hiv.rgc.rr
  hiv.dual.rr <- dat$param$hiv.dual.rr

  # Data
  al <- dat$temp$al


  ## NOTE: In the following section, always force dal into a matrix to avoid
  ##       breaking function when HIV transmission is very low (0 or 1 at-
  ##       risk acts). Possibly the original source of some errors during
  ##       calibration, where extreme parameter sets may have resulted in
  ##       low/no HIV.
  dal <- matrix(
    al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ],
    ncol = ncol(al),
    dimnames = list(NULL, colnames(al))
  )

  if (nrow(dal) == 0) {
    return(dat)
  }

  if (nrow(dal) == 1) {
    dal <- matrix(dal, nrow = 1, dimnames = list(NULL, colnames(dal)))
  } else {
    dal <- dal[sample(1:nrow(dal)), ]
  }

  ncols <- dim(dal)[2]

  ## Reorder by role: ins on the left, rec on the right
  disc.ip <- dal[dal[, "ins"] == 1, ]
  if (is.vector(disc.ip)) {
    disc.ip <- matrix(disc.ip, nrow = 1, dimnames = list(NULL, names(disc.ip)))
  }

  disc.rp <- dal[dal[, "ins"] == 0, c(2:1, 3:ncols)]
  if (is.vector(disc.rp)) {
    disc.rp <- matrix(disc.rp, nrow = 1, dimnames = list(NULL, names(disc.rp)))
  }

  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("ins", "rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]
  ip.rGC <- rGC[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- pmin(0.99, URAI.prob * 2.45^(ip.vl - 4.5))

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1 - ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  condom.rr <- rep(NA, nrow(disc.ip))
  races <- sort(unique(race[disc.ip[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.ip[, 1]] == i))
    condom.rr[not.UAI.race] <- 1 - (cond.eff - cond.fail[i])
  }
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr[not.UAI])

  # PrEP, by adherence class
  ip.on.prep <- which(ip.prep == 1)
  ip.tlo[ip.on.prep] <- ip.tlo[ip.on.prep] + log(prep.hr[ip.prepcl[ip.on.prep]])

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% 1:2)
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.rGC <- which(ip.rGC == 1)
  ip.tlo[is.rGC] <- ip.tlo[is.rGC] + log(hiv.rgc.rr)

  # Race-specific scalar for calibration
  races <- race[disc.ip[, 2]]
  ip.tlo <- ip.tlo + log(trans.scale[races])

  # Convert back to probability
  ip.tprob <- plogis(ip.tlo)
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)


  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]
  rp.uGC <- uGC[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- pmin(0.99, UIAI.prob * 2.45^(rp.vl - 4.5))

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1 - rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  condom.rr <- rep(NA, nrow(disc.rp))
  races <- sort(unique(race[disc.rp[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.rp[, 1]] == i))
    condom.rr[not.UAI.race] <- 1 - (cond.eff - cond.fail[i])
  }
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr[not.UAI])

  # PrEP, by adherence class
  rp.on.prep <- which(rp.prep == 1)
  rp.tlo[rp.on.prep] <- rp.tlo[rp.on.prep] + log(prep.hr[rp.prepcl[rp.on.prep]])

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% 1:2)
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.uGC <- which(rp.uGC == 1)
  rp.tlo[is.uGC] <- rp.tlo[is.uGC] + log(hiv.ugc.rr)

  # Race-specific scalar for calibration
  races <- race[disc.rp[, 1]]
  rp.tlo <- rp.tlo + log(trans.scale[races])

  # Convert back to probability
  rp.tprob <- plogis(rp.tlo)
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)


  # Transmission --------------------------------------------------------

  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------

  infected <- NULL
  if (sum(trans.ip, trans.rp) > 0) {
    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])

    # Attributes of newly infected
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- 1
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
    dat$attr$cuml.time.on.tx[infected] <- 0
    dat$attr$cuml.time.off.tx[infected] <- 0

    # Attributes of transmitter
    transmitter <- c(disc.ip[trans.ip == 1, 1],
                     disc.rp[trans.rp == 1, 2])
    tab.trans <- table(transmitter)
    uni.trans <- as.numeric(names(tab.trans))
    dat$attr$count.trans[uni.trans] <- dat$attr$count.trans[uni.trans] +
                                       as.numeric(tab.trans)
  }

  # Summary Output
  dat$epi$incid[at] <- length(infected)

  ## race/ethnicity-specific
  dat$epi$incid.B[at] <- sum(dat$attr$race[infected] == 1, na.rm = TRUE)
  dat$epi$incid.H[at] <- sum(dat$attr$race[infected] == 2, na.rm = TRUE)
  dat$epi$incid.O[at] <- sum(dat$attr$race[infected] == 3, na.rm = TRUE)
  dat$epi$incid.W[at] <- sum(dat$attr$race[infected] == 4, na.rm = TRUE)

  ## age group-specific
  dat$epi$incid.age1[at] <- sum(dat$attr$age.grp[infected] == 1, na.rm = TRUE)
  dat$epi$incid.age2[at] <- sum(dat$attr$age.grp[infected] == 2, na.rm = TRUE)
  dat$epi$incid.age3[at] <- sum(dat$attr$age.grp[infected] == 3, na.rm = TRUE)
  dat$epi$incid.age4[at] <- sum(dat$attr$age.grp[infected] == 4, na.rm = TRUE)
  dat$epi$incid.age5[at] <- sum(dat$attr$age.grp[infected] == 5, na.rm = TRUE)

  return(dat)
}



hughes_tp <- function(vls, susmales, susages, suscircs, prop.male, fmat = FALSE) {

  suscircs[is.na(suscircs)] <- 0

  sus.hsv2 <- 0.59*prop.male + 0.86*(1 - prop.male)
  sus.gud <- 0.039*prop.male + 0.053*(1 - prop.male)
  sus.tvagin <- 0.068*prop.male + 0.12*(1 - prop.male)
  sus.cerv <- 0.066*(1 - prop.male)

  interc <- -8.3067
  coef.vl <- 1.062566
  coef.male <- 0.6430989
  coef.age <- -0.0403451
  coef.hsv2 <- 0.7625081
  coef.circ <- -0.6377294
  coef.gud <- 0.9749536
  coef.vagin <- 0.9435334
  coef.cerv <- 1.288279

  tp.full <- exp(interc + coef.vl*(vls - 4) +
                   coef.male*susmales + coef.age*(susages - 35) +
                   coef.hsv2*sus.hsv2 + coef.circ*susmales*suscircs +
                   coef.gud*sus.gud + coef.vagin*sus.tvagin +
                   coef.cerv*sus.cerv)

  if (fmat == TRUE) {
    tp.full <- data.frame(tp.full, vls, susmales, susages, suscircs)
  }

  return(tp.full)
}

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}

nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}
