
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
browser()
  # Attributes
  race <- dat$attr$race
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status

  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # Parameters
  rcomp.prob <- dat$param$rcomp.prob
  rcomp.adh.groups <- dat$param$rcomp.adh.groups

  # Condom Use Models
  cond.mc.mod <- param$cond.mc.mod
  cond.oo.mod <- param$cond.oo.mod

  # Temp edgelist
  el <- dat$temp$el

  ## Main/casual partnerships ##
  el.mc <- el[el[, "ptype"] != 3, ]

  race.combo <- race[el.mc[, 1]] + race[el.mc[, 2]]
  comb.age <- age[el.mc[, 1]] + age[el.mc[, 2]]
  hiv.concord.pos <- rep(0, nrow(el.mc))
  cp <- which(diag.status[el.mc[, 1]] == 1 & diag.status[el.mc[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  x <- data.frame(ptype = el.mc[, "ptype"],
                  duration = el.mc[, "durations"],
                  race.combo = race.combo,
                  comb.age = comb.age,
                  hiv.concord.pos = hiv.concord.pos,
                  city = dat$param$netstats$demog$city)
  cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
  el.mc <- cbind(el.mc, cond.prob)

  ## One-off partnerships ##
  el.oo <- el[el[, "ptype"] == 3, ]

  # TODO: consolidate these calcs in el
  race.combo <- race[el.oo[, 1]] + race[el.oo[, 2]]
  comb.age <- age[el.oo[, 1]] + age[el.oo[, 2]]
  hiv.concord.pos <- rep(0, nrow(el.oo))
  cp <- which(diag.status[el.oo[, 1]] == 1 & diag.status[el.oo[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  x <- data.frame(race.combo = race.combo,
                  comb.age = comb.age,
                  hiv.concord.pos = hiv.concord.pos,
                  city = dat$param$netstats$demog$city)
  cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
  el.oo <- cbind(el.oo, cond.prob)


  ## Bind el together
  el <- rbind(el.mc, el.oo)

  # Acts
  ai.vec <- el[, "ai"]
  pid <- rep(1:length(ai.vec), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)
  cond.prob <- rep(el[, "cond.prob"], ai.vec)

  # PrEP Status (risk compensation)
  # if (rcomp.prob > 0) {
  #   idsRC <- which((prepStat[elt[, 1]] == 1 & prepClass[elt[, 1]] %in% rcomp.adh.groups) |
  #                    (prepStat[elt[, 2]] == 1 & prepClass[elt[, 2]] %in% rcomp.adh.groups))
  #   uai.prob[idsRC] <- 1 - (1 - uai.prob[idsRC]) * (1 - rcomp.prob)
  # }

  uai <- rbinom(length(cond.prob), 1, 1 - cond.prob)

  # Act list construction
  al <- cbind(p1, p2, ptype, uai, pid)
  dat$temp$al <- al

  return(dat)
}
