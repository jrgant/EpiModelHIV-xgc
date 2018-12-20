
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {

  # Attributes
  race <- dat$attr$race
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # Parameters
  rcomp.prob <- dat$param$rcomp.prob
  rcomp.adh.groups <- dat$param$rcomp.adh.groups

  el <- dat$temp$el

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Parameters
    cond.rr <- dat$param$cond.rr

    if (type == "main") {
      cond.prob <- dat$param$cond.main.prob
      cond.always <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.prob <- dat$param$cond.pers.prob
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "inst") {
      cond.prob <- dat$param$cond.inst.prob
      cond.always <- dat$attr$cond.always.inst
      ptype <- 3
    }

    elt <- el[el[, "ptype"] == ptype, ]

    ## Process ##

    # Base condom probs
    race.p1 <- race[elt[, 1]]
    race.p2 <- race[elt[, 2]]
    num.B <- (race.p1 == "B") + (race.p2 == "B")
    cond.prob <- (num.B == 2) * (cond.prob[1] * cond.rr[1]) +
                 (num.B == 1) * (cond.prob[2] * cond.rr[2]) +
                 (num.B == 0) * (cond.prob[3] * cond.rr[3])
    uai.prob <- 1 - cond.prob

    # PrEP Status (risk compensation)
    if (rcomp.prob > 0) {
      idsRC <- which((prepStat[elt[, 1]] == 1 & prepClass[elt[, 1]] %in% rcomp.adh.groups) |
                     (prepStat[elt[, 2]] == 1 & prepClass[elt[, 2]] %in% rcomp.adh.groups))
      uai.prob[idsRC] <- 1 - (1 - uai.prob[idsRC]) * (1 - rcomp.prob)
    }

    ai.vec <- elt[, "ai"]
    p1 <- rep(elt[, "p1"], ai.vec)
    p2 <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(p1), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(p1, p2, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(p1, p2, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al


  return(dat)
}
