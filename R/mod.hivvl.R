
#' @title Viral Load Module
#'
#' @description Module function for updating HIV viral load.
#'
#' @inheritParams aging_msm
#'
#' @details
#' HIV viral load varies over time as a function of time since infection and ART
#' history. In the absence of ART, VL rises during the acute rising stage and
#' falls during the acute falling stage, until it reaches a set-point value in
#' chronic stage infection. VL again rises during AIDS stage disease until the
#' point of death.
#'
#' For persons who have ever initated treatment (\code{tt.traj} is \code{3} or
#' \code{4}), VL changes depending on current ART use in that time step.
#' Current use is associated with a reduction in VL, with the rates of decline
#' and nadirs dependent on partial or full suppression levels. Current
#' non-adherence is associated with an equal level of increase to VL. All persons
#' who have reached AIDS, regardless of how they arrived, have a similar rate of
#' VL increase.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{vl} attribute.
#'
#' @keywords module msm
#'
#' @export
#'
hivvl_msm <- function(dat, at) {

  # Attributes
  time.inf <- at - dat$attr$inf.time
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  status <- dat$attr$status
  tt.traj <- dat$attr$tt.traj
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  tx.status <- dat$attr$tx.status

  # Parameters
  acute.rise.int <- dat$param$vl.acute.rise.int
  acute.peak <- dat$param$vl.acute.peak
  acute.fall.int <- dat$param$vl.acute.fall.int
  vl.set.point <- dat$param$vl.set.point
  aids.onset <- dat$param$vl.aids.onset
  aids.int <- dat$param$vl.aids.int
  vl.fatal <- dat$param$vl.fatal
  vl.full.supp <- dat$param$vl.full.supp
  vl.tx.down.slope <- dat$param$vl.tx.down.slope
  vl.part.supp <- dat$param$vl.part.supp
  vl.tx.up.slope <- dat$param$vl.tx.up.slope
  vl.aids.slope <- (vl.fatal - vl.set.point) / aids.int

  ## Process

  # 1. TX-naive
  idsElig1 <- which(status == 1 & cum.time.on.tx == 0)
  time.inf1 <- time.inf[idsElig1]
  new.vl <- rep(NA, length(idsElig1))

  # Acute rising
  idsElig1.AR <- which(stage[idsElig1] == 1)
  new.vl[idsElig1.AR] <- acute.peak * time.inf1[idsElig1.AR] / acute.rise.int

  # Acute falling
  idsElig1.AF <- which(stage[idsElig1] == 2)
  new.vl[idsElig1.AF] <- ((vl.set.point - acute.peak) *
                            (time.inf1[idsElig1.AF] - acute.rise.int) / acute.fall.int + acute.peak)

  # Chronic
  idsElig1.C <- which(stage[idsElig1] == 3)
  new.vl[idsElig1.C] <- vl.set.point

  # AIDS
  idsElig1.A <- which(stage[idsElig1] == 4)
  new.vl[idsElig1.A] <- vl.set.point + (time.inf1[idsElig1.A] - aids.onset) * vl.aids.slope

  vl[idsElig1] <- new.vl

  # 2. On tx, tt.traj=full/dur, not AIDS
  target <- which(tx.status == 1 & tt.traj %in% 2:3 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - vl.tx.down.slope, vl.full.supp)
  vl[target] <- new.vl

  # 3. On tx, tt.traj=part, not AIDS
  target <- which(tx.status == 1 & tt.traj == 1 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - vl.tx.down.slope, vl.part.supp)
  vl[target] <- new.vl

  # 4. Off tx, not naive, tt.traj=full/dur, not AIDS
  target <- which(tx.status == 0 & tt.traj %in% 2:3 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + vl.tx.up.slope, vl.set.point)
  vl[target] <- new.vl

  # 5. Off tx, not naive, tt.traj=part, not AIDS
  target <- which(tx.status == 0 & tt.traj == 1 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + vl.tx.up.slope, vl.set.point)
  vl[target] <- new.vl

  # 6. On tx, tt.traj=full/dur, AIDS
  # NA

  # 7. On tx, tt.traj=part, AIDS (check this group reduces VL to set point)
  target <- which(tx.status == 1 &
                  tt.traj == 1 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vl.aids.slope
  vl[target] <- new.vl

  # 8. Off tx, tt.traj=full/dur and AIDS
  target <- which(tx.status == 0 & tt.traj %in% 2:3 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vl.aids.slope
  vl[target] <- new.vl

  # 9. Off tx, tt.traj=part, and AIDS (check this group increases VL to right level)
  target <- which(tx.status == 0 & tt.traj == 1 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vl.aids.slope
  vl[target] <- new.vl


  ## Output
  dat$attr$vl <- vl

  idsSupp <- which(vl <= log10(200))
  idsUsupp <- which(vl > log10(200))
  dat$attr$vl.last.usupp[idsUsupp] <- at
  dat$attr$vl.last.supp[idsSupp] <- at

  if (dat$control$save.clin.hist == TRUE) {
    dat <- save_clin_hist(dat, at)
  }

  return(dat)
}

save_clin_hist <- function(dat, at) {

  if (is.null(dat$temp$clin.hist)) {
    m <- list()
    for (i in 1:3) {
      m[[i]] <- array(dim = c(length(dat$attr$active), dat$control$nsteps))
    }
  } else {
    m <- dat$temp$clin.hist
  }
  m[[1]][, at] <- dat$attr$vl
  m[[2]][, at] <- dat$attr$stage
  m[[3]][, at] <- dat$attr$tx.status

  dat$temp$clin.hist <- m
  return(dat)
}


#' @export
#' @rdname hivvl_msm
vl_het <- function(dat, at) {

  ## Common variables
  status <- dat$attr$status
  infTime <- dat$attr$infTime


  # Assign base VL ----------------------------------------------------------
  if (is.null(dat$attr$vlLevel)) {
    dat$attr$vlLevel <- rep(NA, length(status))
    dat$attr$vlSlope <- rep(NA, length(status))
  }
  vlLevel <- dat$attr$vlLevel

  idsEligAsn <- which(status == 1 & is.na(vlLevel))
  if (length(idsEligAsn) > 0) {
    vlLevel[idsEligAsn] <- expected_vl(male = dat$attr$male[idsEligAsn],
                                       age = dat$attr$age[idsEligAsn],
                                       ageInf = dat$attr$ageInf[idsEligAsn],
                                       param = dat$param)
  }


  # Update natural VL -------------------------------------------------------
  txStartTime <- dat$attr$txStartTime
  idsEligUpd <- which(status == 1 &
                      infTime < at & is.na(txStartTime))

  if (length(idsEligUpd) > 0) {
    vlLevel[idsEligUpd] <- expected_vl(male = dat$attr$male[idsEligUpd],
                                       age = dat$attr$age[idsEligUpd],
                                       ageInf = dat$attr$ageInf[idsEligUpd],
                                       param = dat$param)
  }

  # VL decline with ART -----------------------------------------------------
  txStat <- dat$attr$txStat
  idsEligTx <- which(status == 1 & infTime < at & txStat == 1)
  if (length(idsEligTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time
    tx.vlsupp.level <- dat$param$tx.vlsupp.level

    vlSlope <- dat$attr$vlSlope
    needSlope <- intersect(idsEligTx, which(is.na(vlSlope)))

    vl.slope <- vlSlope
    if (length(needSlope) > 0) {
      vl.diff <- pmin(tx.vlsupp.level - vlLevel[needSlope], 0)
      vl.slope[needSlope] <- vl.diff / tx.vlsupp.time
      dat$attr$vlSlope[needSlope] <- vl.slope[needSlope]
    }

    vlLevel[idsEligTx] <- pmax(vlLevel[idsEligTx] + vl.slope[idsEligTx], tx.vlsupp.level)
  }


  # VL rebound post ART -----------------------------------------------------
  idsEligNoTx <- which(status == 1 &
                       txStat == 0 & !is.na(txStartTime))
  if (length(idsEligNoTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time

    expVl <- expected_vl(male = dat$attr$male[idsEligNoTx],
                         age = dat$attr$age[idsEligNoTx],
                         ageInf = dat$attr$ageInf[idsEligNoTx],
                         param = dat$param)

    vl.slope <- dat$attr$vlSlope

    vlLevel[idsEligNoTx] <- pmin(vlLevel[idsEligNoTx] - vl.slope[idsEligNoTx], expVl)
  }

  dat$attr$vlLevel <- vlLevel

  return(dat)
}


expected_vl <- function(male, age, ageInf, param) {

  timeInf <- (age - ageInf) * (365 / param$time.unit)

  slope1 <- param$vl.acute.peak / param$vl.acute.topeak
  slope2 <- (param$vl.setpoint - param$vl.acute.peak) /
    (param$vl.acute.toset - param$vl.acute.topeak)

  sl3denom <- expected_cd4(method = "timeto",
                           cd4Count1 = 200, cd4Count2 = 25,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)
  slope3 <- (param$vl.aidsmax - param$vl.setpoint) / sl3denom

  setptTime <- param$vl.acute.topeak + param$vl.acute.toset
  aidsTime <- expected_cd4(method = "timeto", cd4Count1 = 200,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)

  gp <- 1 * (timeInf <= param$vl.acute.topeak) +
    2 * (timeInf > param$vl.acute.topeak & timeInf <= setptTime) +
    3 * (timeInf > setptTime & timeInf <= aidsTime) +
    4 * (timeInf > aidsTime)

  vlLevel <- rep(NA, length(timeInf))
  vlLevel[gp == 1] <- timeInf[gp == 1] * slope1
  vlLevel[gp == 2] <- pmax(param$vl.setpoint,
                           param$vl.acute.peak +
                             (timeInf[gp == 2] - param$vl.acute.topeak) * slope2)
  vlLevel[gp == 3] <- param$vl.setpoint
  vlLevel[gp == 4] <- pmin(param$vl.aidsmax,
                           param$vl.setpoint +
                             (timeInf[gp == 4] - aidsTime[gp == 4]) * slope3[gp == 4])

  return(vlLevel)
}
