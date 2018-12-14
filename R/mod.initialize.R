
# MSM -----------------------------------------------------------------

#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module msm
#'
initialize_msm <- function(x, param, init, control, s) {

  # Master data list
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  assignInNamespace("InitErgmConstraint..attributes",
                    EpiModelHIV::InitErgmConstraint..attributes,
                    ns = "ergm", envir = as.environment("package:ergm"))

  ## Network simulation ##
  nw <- list()
  for (i in 1:3) {
    nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }

  ## Build initial edgelists
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    attributes(dat$el[[i]])$vnames <- NULL
    p <- stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                     x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- tergmLite::ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p


  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }


  ## Nodal attributes ##

  # Degree terms
  dat$attr$deg.pers <- get.vertex.attribute(x[[1]]$fit$network, "deg.pers")
  dat$attr$deg.main <- get.vertex.attribute(x[[2]]$fit$network, "deg.main")


  # Race
  dat$attr$race <- get.vertex.attribute(nw[[1]], "race")
  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")

  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$sqrt.age <- get.vertex.attribute(nw[[1]], "sqrt.age")
  dat$attr$age <- dat$attr$sqrt.age^2

  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[3]], "riskg")

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(num, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers <- uai.always[, 1]
  dat$attr$cond.always.inst <- uai.always[, 2]

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(NA, num)
  circ[ids.B] <- sample(apportion_lr(num.B, 0:1, 1 - param$circ.prob[1]))
  circ[ids.W] <- sample(apportion_lr(num.W, 0:1, 1 - param$circ.prob[2]))
  dat$attr$circ <- circ

  ## PrEP Attributes ##
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

  # Role class
  role.class <- get.vertex.attribute(nw[[1]], "role.class")
  dat$attr$role.class <- role.class

  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[role.class == "I"]  <- 1
  ins.quot[role.class == "R"]  <- 0
  ins.quot[role.class == "V"]  <- runif(sum(role.class == "V"))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msm(dat)

  ## GC/CT status
  idsUreth <- which(role.class %in% c("I", "V"))
  idsRect <- which(role.class %in% c("R", "V"))

  uGC <- rGC <- rep(0, num)
  uCT <- rCT <- rep(0, num)

  # Initialize GC infection at both sites
  idsUGC <- sample(idsUreth, size = round(init$prev.ugc * num), FALSE)
  uGC[idsUGC] <- 1

  idsRGC <- sample(setdiff(idsRect, idsUGC), size = round(init$prev.rgc * num), FALSE)
  rGC[idsRGC] <- 1

  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC

  dat$attr$rGC.sympt <- dat$attr$uGC.sympt <- rep(NA, num)
  dat$attr$rGC.sympt[rGC == 1] <- rbinom(sum(rGC == 1), 1, dat$param$rgc.sympt.prob)
  dat$attr$uGC.sympt[uGC == 1] <- rbinom(sum(uGC == 1), 1, dat$param$ugc.sympt.prob)

  dat$attr$rGC.infTime <- dat$attr$uGC.infTime <- rep(NA, length(dat$attr$active))
  dat$attr$rGC.infTime[rGC == 1] <- 1
  dat$attr$uGC.infTime[uGC == 1] <- 1

  dat$attr$rGC.timesInf <- rep(0, num)
  dat$attr$rGC.timesInf[rGC == 1] <- 1
  dat$attr$uGC.timesInf <- rep(0, num)
  dat$attr$uGC.timesInf[uGC == 1] <- 1

  dat$attr$rGC.tx <- dat$attr$uGC.tx <- rep(NA, num)
  dat$attr$rGC.tx.prep <- dat$attr$uGC.tx.prep <- rep(NA, num)

  # Initialize CT infection at both sites
  idsUCT <- sample(idsUreth, size = round(init$prev.uct * num), FALSE)
  uCT[idsUCT] <- 1

  idsRCT <- sample(setdiff(idsRect, idsUCT), size = round(init$prev.rct * num), FALSE)
  rCT[idsRCT] <- 1

  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rCT.sympt <- dat$attr$uCT.sympt <- rep(NA, num)
  dat$attr$rCT.sympt[rCT == 1] <- rbinom(sum(rCT == 1), 1, dat$param$rct.sympt.prob)
  dat$attr$uCT.sympt[uCT == 1] <- rbinom(sum(uCT == 1), 1, dat$param$uct.sympt.prob)

  dat$attr$rCT.infTime <- dat$attr$uCT.infTime <- rep(NA, num)
  dat$attr$rCT.infTime[dat$attr$rCT == 1] <- 1
  dat$attr$uCT.infTime[dat$attr$uCT == 1] <- 1

  dat$attr$rCT.timesInf <- rep(0, num)
  dat$attr$rCT.timesInf[rCT == 1] <- 1
  dat$attr$uCT.timesInf <- rep(0, num)
  dat$attr$uCT.timesInf[uCT == 1] <- 1

  dat$attr$rCT.tx <- dat$attr$uCT.tx <- rep(NA, num)
  dat$attr$rCT.tx.prep <- dat$attr$uCT.tx.prep <- rep(NA, num)


  # Network statistics
  dat$stats$nwstats <- list()


  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat <- prevalence_msm(dat, at = 1)

  # Save partner list
  plist <- cbind(dat$el[[1]], ptype = 1)
  plist <- rbind(plist, cbind(dat$el[[2]], ptype = 2))
  plist <- cbind(plist, start = 1, stop = NA)
  colnames(plist)[1:2] <- c("p1", "p2")
  dat$temp$plist <- plist

  class(dat) <- "dat"
  return(dat)
}


#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_status_msm <- function(dat) {

  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")
  age <- dat$attr$age
  race <- dat$attr$race

  # Infection Status
  nInfB <- round(dat$init$prev.B * num.B)
  nInfW <- round(dat$init$prev.W * num.W)

  # Age-based infection probability
  probInfCrB <- age[ids.B] * dat$init$init.prev.age.slope.B
  probInfB <- probInfCrB + (nInfB - sum(probInfCrB)) / num.B

  probInfCrW <- age[ids.W] * dat$init$init.prev.age.slope.W
  probInfW <- probInfCrW + (nInfW - sum(probInfCrW)) / num.W

  if (any(probInfB <= 0) | any(probInfW <= 0)) {
    stop("Slope of initial prevalence by age must be sufficiently low to ",
         "avoid non-positive probabilities.", call. = FALSE)
  }

  # Infection status
  status <- rep(0, num)
  while (sum(status[ids.B]) != nInfB) {
    status[ids.B] <- rbinom(num.B, 1, probInfB)
  }
  while (sum(status[ids.W]) != nInfW) {
    status[ids.W] <- rbinom(num.W, 1, probInfW)
  }
  dat$attr$status <- status


  # Treatment trajectory
  tt.traj <- rep(NA, num)

  tt.traj[ids.B] <- sample(apportion_lr(num.B, 1:4,
                                        dat$param$tt.traj.prob[[1]]))
  tt.traj[ids.W] <- sample(apportion_lr(num.W, 1:4,
                                        dat$param$tt.traj.prob[[2]]))
  dat$attr$tt.traj <- tt.traj



  ## Infection-related attributes

  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  last.neg.test <- rep(NA, num)
  tx.status <- rep(NA, num)
  tx.init.time <- rep(NA, num)
  cum.time.on.tx <- rep(NA, num)
  cum.time.off.tx <- rep(NA, num)
  infector <- rep(NA, num)
  inf.role <- rep(NA, num)
  inf.type <- rep(NA, num)
  inf.diag <- rep(NA, num)
  inf.tx <- rep(NA, num)
  inf.stage <- rep(NA, num)

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                  min(dat$init$ages), 0))

  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vlds <- (vlf - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int


  ### Non-treater type: tester and non-tester
  selected <- which(status == 1 & tt.traj %in% c(1, 2))
  max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  tx.status[selected] <- 0
  cum.time.on.tx[selected] <- 0
  cum.time.off.tx[selected] <- time.since.inf

  stage[selected[time.since.inf <= vlar.int]] <- 1
  stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- 2
  stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- 3
  stage[selected[time.since.inf > vldo.int]] <- 4

  stage.time[selected][stage[selected] == 1] <- time.since.inf[stage[selected] == 1]
  stage.time[selected][stage[selected] == 2] <- time.since.inf[stage[selected] == 2] -
                                                   vlar.int
  stage.time[selected][stage[selected] == 3] <- time.since.inf[stage[selected] == 3] -
                                                  vl.acute.int
  stage.time[selected][stage[selected] == 4] <- time.since.inf[stage[selected] == 4] -
                                                  vldo.int

  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
                  (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

  selected <- which(status == 1 & tt.traj == 1)
  diag.status[selected] <- 0

  selected <- which(status == 1 & tt.traj == 2)

  # Time to next test
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$hiv.test.int[1] * (race[selected] == "B") +
                        dat$param$hiv.test.int[2] * (race[selected] == "W")))

  twind.int <- dat$param$test.window.int
  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1


  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

  tx.init.time.B <- twind.int + dat$param$hiv.test.int[1] + 1 / dat$param$tx.init.prob[1]
  tx.init.time.W <- twind.int + dat$param$hiv.test.int[2] + 1 / dat$param$tx.init.prob[2]

  # Stage for Blacks
  prop.time.on.tx.B <- dat$param$tx.reinit.prob[1] /
                       (dat$param$tx.halt.prob[1] + dat$param$tx.reinit.prob[1])
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  numsteps.B <- (dat$param$max.time.off.tx.full.int - tx.init.time.B) /
                (1 - prop.time.on.tx.B)
  offon.B <- rbind(offon.B,
                   cbind(tx.init.time.B + (1 - prop.time.on.tx.B) * 1:numsteps.B,
                         prop.time.on.tx.B * 1:numsteps.B))
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  # Stage for Whites
  prop.time.on.tx.W <- dat$param$tx.reinit.prob[2] /
    (dat$param$tx.halt.prob[2] + dat$param$tx.reinit.prob[2])
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)
  numsteps.W <- (dat$param$max.time.off.tx.full.int - tx.init.time.W) /
    (1 - prop.time.on.tx.W)
  offon.W <- rbind(offon.W,
                   cbind(tx.init.time.W + (1 - prop.time.on.tx.W) * 1:numsteps.W,
                         prop.time.on.tx.W * 1:numsteps.W))
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # Vl for Blacks
  selected <- which(status == 1 & tt.traj == 4 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 4 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # Diagnosis
  selected <- which(status == 1 & tt.traj == 4)
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$hiv.test.int[1] * (race[selected] == "B") +
                        dat$param$hiv.test.int[2] * (race[selected] == "W")))

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]
  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  ### Part adherent type

  # Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

  prop.time.on.tx.B <- dat$param$tx.reinit.prob[1] /
                       (dat$param$tx.halt.prob[1] + dat$param$tx.reinit.prob[1])
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  while (offon.B[nrow(offon.B), 1] / dat$param$max.time.off.tx.part.int +
         offon.B[nrow(offon.B), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B <- rbind(offon.B,
                     offon.B[nrow(offon.B), ] + c(1 - prop.time.on.tx.B,
                                                      prop.time.on.tx.B))
  }
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  prop.time.on.tx.W <- dat$param$tx.reinit.prob[2] /
                       (dat$param$tx.halt.prob[2] + dat$param$tx.reinit.prob[2])
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)

  while (offon.W[nrow(offon.W), 1] / dat$param$max.time.off.tx.part.int +
         offon.W[nrow(offon.W), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W <- rbind(offon.W,
                     offon.W[nrow(offon.W), ] + c(1 - prop.time.on.tx.W,
                                                  prop.time.on.tx.W))
  }
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # VL for Blacks
  selected <- which(status == 1 & tt.traj == 3 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 3 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == 3)
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$hiv.test.int[1] * (race[selected] == "B") +
                        dat$param$hiv.test.int[2] * (race[selected] == "W")))

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(2, 3, 4))
  tslt <- rgeom(length(selected),
                1 / (dat$param$hiv.test.int[1] * (race[selected] == "B") +
                     dat$param$hiv.test.int[1] * (race[selected] == "W")))
  last.neg.test[selected] <- -tslt


  ## Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$tx.status <- tx.status
  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx
  dat$attr$infector <- infector
  dat$attr$inf.role <- inf.role
  dat$attr$inf.type <- inf.type
  dat$attr$inf.diag <- inf.diag
  dat$attr$inf.tx <- inf.tx
  dat$attr$inf.stage <- inf.stage

  return(dat)

}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module msm
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}


#' @export
InitErgmConstraint..attributes <- function(lhs.nw, ...){
  list(
    # free_dyads = {
    # n <- network.size(lhs.nw)
    ## NB: Free dyad RLE matrix is stored in a column-major order for
    ## consistency with R.
    # d <-
    # if(has.loops(lhs.nw)) rep(rep(rle(TRUE),n,scale="run"),n,scale="run")
    # else do.call(c, lapply(seq_len(n), function(i) rep(rle(c(TRUE,FALSE,TRUE)), c(i-1, 1, n-i),scale="run")))

    # if(is.bipartite(lhs.nw)){
    #   n1 <- lhs.nw%n%"bipartite"
    #   n2 <- n - n1
    #
    #   d <- d &
    #     c(rep(rep(rle(c(FALSE)),n,scale="run"),n1,scale="run"),
    #       rep(rep(rle(c(TRUE,FALSE)),c(n1,n2),scale="run"),n2,scale="run"))
    # }

    # if(!is.directed(lhs.nw)){
    #   d <- d &
    #     do.call(c, lapply(seq_len(n), function(i) rep(rle(c(TRUE,FALSE)), c(i-1, n-i+1),scale="run")))
    # }
    #
    # rlebdm(compact.rle(d), n)
    # },
    constrain = character(0),
    dependence = FALSE)
}


# HET -----------------------------------------------------------------


#' @export
#' @rdname initialize_msm
initialize_het <- function(x, param, init, control, s) {

  dat <- list()
  dat$temp <- list()
  nw <- simulate(x$fit, control = control.simulate.ergm(MCMC.burnin = 1e6))

  dat$el <- list()
  dat$el[[1]] <- as.edgelist(nw)
  attributes(dat$el)$vnames <- NULL
  p <- tergmLite::stergm_prep(nw, x$formation, x$coef.diss$dissolution, x$coef.form,
                              x$coef.diss$coef.adj, x$constraints)
  p$model.form$formula <- NULL
  p$model.diss$formula <- NULL
  dat$p <- list()
  dat$p[[1]] <- p

  ## Network Model Parameters
  dat$nwparam <- list(x[-which(names(x) == "fit")])

  ## Simulation Parameters
  dat$param <- param
  dat$param$modes <- 1

  dat$init <- init
  dat$control <- control

  ## Nodal Attributes
  dat$attr <- list()

  dat$attr$male <- get.vertex.attribute(nw, "male")

  n <- network.size(nw)
  dat$attr$active <- rep(1, n)
  dat$attr$entTime <- rep(1, n)

  dat <- initStatus_het(dat)

  age <- rep(NA, n)
  age[dat$attr$male == 0] <- sample(init$ages.feml, sum(dat$attr$male == 0), TRUE)
  age[dat$attr$male == 1] <- sample(init$ages.male, sum(dat$attr$male == 1), TRUE)
  dat$attr$age <- age

  dat <- initInfTime_het(dat)
  dat <- initDx_het(dat)
  dat <- initTx_het(dat)

  # Circumcision
  male <- dat$attr$male
  nMales <- sum(male == 1)
  age <- dat$attr$age

  circStat <- circTime <- rep(NA, n)

  circStat[male == 1] <- rbinom(nMales, 1, dat$param$circ.prob.birth)

  isCirc <- which(circStat == 1)
  circTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

  dat$attr$circStat <- circStat
  dat$attr$circTime <- circTime


  ## Stats List
  dat$stats <- list()

  ## Final steps
  dat$epi <- list()
  dat <- prevalence_het(dat, at = 1)

}


#' @title Reinitialization Module
#'
#' @description This function reinitializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_het}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_het}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_het}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @keywords module het
#'
#' @export
#'
reinit_het <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi",
                       "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"

  return(dat)
}


initStatus_het <- function(dat) {

  ## Variables
  i.prev.male <- dat$init$i.prev.male
  i.prev.feml <- dat$init$i.prev.feml

  male <- dat$attr$male
  idsMale <- which(male == 1)
  idsFeml <- which(male == 0)
  nMale <- length(idsMale)
  nFeml <- length(idsFeml)
  n <- nMale + nFeml

  ## Process
  status <- rep(0, n)
  status[sample(idsMale, round(i.prev.male * nMale))] <- 1
  status[sample(idsFeml, round(i.prev.feml * nFeml))] <- 1

  dat$attr$status <- status

  return(dat)
}


initInfTime_het <- function(dat) {

  status <- dat$attr$status
  n <- length(status)

  infecteds <- which(status == 1)
  infTime <- rep(NA, n)

  inf.time.dist <- dat$init$inf.time.dist

  if (inf.time.dist == "allacute") {
    max.inf.time <- dat$param$vl.acute.topeak + dat$param$vl.acute.toset
    infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
  } else {
    max.inf.time <- dat$init$max.inf.time / dat$param$time.unit
    if (inf.time.dist == "geometric") {
      total.d.rate <- 1/max.inf.time
      infTime[infecteds] <- -rgeom(length(infecteds), total.d.rate)
    }
    if (inf.time.dist == "uniform") {
      infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
    }
  }

  ## Enforce that time infected < age
  infTime[infecteds] <- pmax(infTime[infecteds],
                             1 - dat$attr$age[infecteds] * (365 / dat$param$time.unit))

  dat$attr$infTime <- infTime

  timeInf <- 1 - infTime
  dat$attr$ageInf <- pmax(0, dat$attr$age - round(timeInf) * (dat$param$time.unit / 365))

  stopifnot(all(dat$attr$ageInf[infecteds] <= dat$attr$age[infecteds]),
            all(dat$attr$ageInf[infecteds] >= 0))

  return(dat)
}


initDx_het <- function(dat) {

  n <- sum(dat$attr$active == 1)
  status <- dat$attr$status

  dxStat <- rep(NA, n)
  dxStat[status == 1] <- 0

  dxTime <- rep(NA, n)

  dat$attr$dxStat <- dxStat
  dat$attr$dxTime <- dxTime

  return(dat)
}


initTx_het <- function(dat) {

  ## Variables
  status <- dat$attr$status
  n <- sum(dat$attr$active == 1)
  nInf <- sum(status == 1)

  tx.init.cd4.mean <- dat$param$tx.init.cd4.mean
  tx.init.cd4.sd <- dat$param$tx.init.cd4.sd
  tx.elig.cd4 <- dat$param$tx.elig.cd4


  ## Process
  dat$attr$txStat <- rep(NA, n)
  dat$attr$txStartTime <- rep(NA, n)
  dat$attr$txStops <- rep(NA, n)
  dat$attr$txTimeOn <- rep(NA, n)
  dat$attr$txTimeOff <- rep(NA, n)

  txCD4min <- rep(NA, n)
  txCD4min[status == 1] <- pmin(rnbinom(nInf,
                                        size = nbsdtosize(tx.init.cd4.mean,
                                                          tx.init.cd4.sd),
                                        mu = tx.init.cd4.mean), tx.elig.cd4)
  dat$attr$txCD4min <- txCD4min
  dat$attr$txCD4start <- rep(NA, n)
  dat$attr$txType <- rep(NA, n)

  return(dat)
}

