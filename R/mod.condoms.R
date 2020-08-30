
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
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat

  # Condom Use Models (anal sex acts only)
  cond.mc.mod <- dat$param$epistats$cp.mod
  cond.oo.mod <- dat$param$epistats$cp.mod

  cond.scale <- dat$param$cond.scale

  # Temp edgelist
  el <- dat$temp$el
  el.mc <- el[el[, "ptype"] != 3, ]

  # Get ego/partner race combination
  race.combo <- rep(NA, nrow(el.mc))

  for (i in 1:length(race.combo)) {
    race.combo[i] <- paste0(
      sort(c(race[el.mc[i, 1]], race[el.mc[i, 2]])),
      collapse = ""
    )
  }


  # Get ego/partner age group combination
  age.el1 <- rep(NA, nrow(el.mc))
  age.el2 <- rep(NA, nrow(el.mc))
  age.combo <- rep(NA, nrow(el.mc))
  age.breaks <- dat$param$netstats$demog$age.breaks

  age.el1 <- cut(
    age[el.mc[, 1]],
    age.breaks,
    right = FALSE,
    labels = FALSE
  )

  age.el2 <- cut(
    age[el.mc[, 2]],
    age.breaks,
    right = FALSE,
    labels = FALSE
  )

  if (!(all.equal(length(age.combo), length(age.el1), length(age.el2)))) {
    stop("age.combo, age.el1, and age.el2 must all be the same length.")
  }

  for (i in 1:length(age.combo)) {
    age.combo[i] <- paste0(
      sort(c(age.el1[i], age.el2[i])),
      collapse = ""
    )
  }

  ## HIV concordance within partnerships (all ptypes)
  hiv.concord <- rep(0, nrow(el))
  bothneg <- which(diag.status[el[, 1]] == 0 & diag.status[el[, 2]] == 0)
  bothpos <- which(diag.status[el[, 1]] == 1 & diag.status[el[, 2]] == 1)
  discord <- which(diag.status[el[, 1]] != diag.status[el[, 2]])

  hiv.concord[bothneg] <- 1
  hiv.concord[bothpos] <- 2
  hiv.concord[discord] <- 3

  any.prep <- as.numeric((prepStat[el[, 1]] + prepStat[el[, 2]]) > 0)

  ## Main/casual partnerships ##
  mc.parts <- which(el[, "ptype"] != 3)
  el.mc <- el[mc.parts, ]

  x <- data.frame(
    ptype = el.mc[, "ptype"],
    duration = el.mc[, "durations"],
    race.combo = race.combo[mc.parts],
    age.i = age[el.mc[, 1]],
    age.j = age[el.mc[, 2]],
    hiv.concord = hiv.concord[mc.parts],
    prep = any.prep[mc.parts]
  )

  cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
  el.mc <- cbind(el.mc, cond.prob)

  ## One-off partnerships ##
  oo.parts <- which(el[, "ptype"] == 3)
  el.oo <- el[oo.parts, ]

  x <- data.frame(
    race.combo = race.combo[oo.parts],
    age.i = el.oo[, 1],
    age.j = el.oo[, 2],
    hiv.concord = hiv.concord[oo.parts],
    prep = any.prep[oo.parts]
  )

  cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
  el.oo <- cbind(el.oo, cond.prob)

  ## Bind el together
  el <- rbind(el.mc, el.oo)

  # anal sex acts
  ai.vec <- el[, "ai"]
  pid <- rep(1:length(ai.vec), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)
  cond.prob <- rep(el[, "cond.prob"], ai.vec)

  cond.prob <- cond.prob * cond.scale

  ## Draw for whether anal sex act is unprotected.
  uai <- rbinom(length(cond.prob), 1, 1 - cond.prob)

  # Anal act list construction
  al <- cbind(p1, p2, ptype, uai, pid)
  dat$temp$al <- al

  # Oral act list construction
  pid <- p1 <- p2 <- ptype <- NULL
  oi.vec <- el[, "oi"]
  pid <- rep(seq_len(length(oi.vec)), oi.vec)
  p1 <- rep(el[, "p1"], oi.vec)
  p2 <- rep(el[, "p2"], oi.vec)
  ptype <- rep(el[, "ptype"], oi.vec)

  ol <- cbind(p1, p2, ptype, pid)
  dat$temp$ol <- ol

  return(dat)
}
