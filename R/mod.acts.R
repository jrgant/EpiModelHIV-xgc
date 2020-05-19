
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and the expected values within black-black, black-white,
#' and white-white combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @keywords module msm
#' @export
#'
acts_msm <- function(dat, at) {

  # Attributes
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  age <- dat$attr$age
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  uid <- dat$attr$uid

  plist <- dat$temp$plist

  # Parameters
  ai.acts.mod <- dat$param$epistats$ai.acts.mod
  oi.acts.mod <- dat$param$epistats$oi.acts.mod
  acts.aids.vl <- dat$param$acts.aids.vl
  ai.acts.scale <- dat$param$ai.acts.scale
  oi.acts.scale <- dat$param$oi.acts.scale

  # Construct edgelist
  el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])
  ptype  <- rep(1:3, times = c(nrow(dat$el[[1]]),
                               nrow(dat$el[[2]]),
                               nrow(dat$el[[3]])))
  st1 <- status[el[, 1]]
  st2 <- status[el[, 2]]

  el <- cbind(el, st1, st2, ptype)
  colnames(el) <- c("p1", "p2", "st1", "st2", "ptype")

  # Subset to main/casual
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

  # Current partnership durations
  pid_plist <- plist[, 1] * 1e7 + plist[, 2]
  pid_el <- uid[el.mc[, 1]] * 1e7 + uid[el.mc[, 2]]
  matches <- match(pid_el, pid_plist)
  durations <- (at - plist[, "start"])[matches]

  # HIV-positive concordant
  hiv.concord.pos <- rep(0, nrow(el.mc))
  cp <- which(diag.status[el.mc[, 1]] == 1 & diag.status[el.mc[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  # Model predictions
  x <- data.frame(
    ptype = el.mc[, "ptype"],
    durat_wks = durations,
    race.combo = race.combo,
    age.combo = age.combo,
    hiv.concord.pos = hiv.concord.pos
  )

  # Predict anal act rates
  ai.rates <- unname(
    predict(ai.acts.mod, newdata = x, type = "response")
  ) / 52

  ai.rates <- ai.rates * ai.acts.scale
  ai <- rpois(length(ai.rates), ai.rates)
  el.mc <- cbind(el.mc, durations, ai)

  # Predict oral act rates
  oi.rates <- unname(
    predict(oi.acts.mod, newdata = x, type = "response")
  ) / 52

  oi.rates <- oi.rates * oi.acts.scale
  oi <- rpois(length(oi.rates), oi.rates)
  el.mc <- cbind(el.mc, oi)

  # Add one-time partnerships
  el.oo <- el[el[, "ptype"] == 3, ]
  ai <- oi <- durations <- rep(1, nrow(el.oo))
  el.oo <- cbind(el.oo, durations, ai, oi)

  # Bind el back together
  el <- rbind(el.mc, el.oo)

  # For AIDS cases with VL above acts.aids.vl, reduce their their acts to 0
  p1HIV <- which(el[, "st1"] == 1)
  p1AIDS <- stage[el[p1HIV, "p1"]] == 4 & vl[el[p1HIV, "p1"]] >= acts.aids.vl
  el[p1HIV[p1AIDS == TRUE], "ai"] <- 0

  p2HIV <- which(el[, "st2"] == 1)
  p2AIDS <- stage[el[p2HIV, "p2"]] == 4 & vl[el[p2HIV, "p2"]] >= acts.aids.vl
  el[p2HIV[p2AIDS == TRUE], "ai"] <- 0

  # Flip order of discordant edges
  disc <- abs(el[, "st1"] - el[, "st2"]) == 1
  disc.st2pos <- which(disc == TRUE & el[, "st2"] == 1)
  el[disc.st2pos, 1:4] <- el[disc.st2pos, c(2, 1, 4, 3)]

  # Remove inactive edges from el
  el <- el[-which(el[, "ai"] == 0), ]

  # Save out
  dat$temp$el <- el

  return(dat)
}
