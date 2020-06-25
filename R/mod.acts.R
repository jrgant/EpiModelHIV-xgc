
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
  age.i <- rep(NA, nrow(el.mc))
  age.j <- rep(NA, nrow(el.mc))

  age.i <- age[el.mc[, 1]]
  age.j <- age[el.mc[, 2]]
  abs_sqrt_agediff <- abs(sqrt(age.i) - sqrt(age.j))

  if (!length(age.i) == length(age.j)) {
    stop("age.i, and age.j must all be the same length.")
  }

  # Current partnership durations
  pid_plist <- plist[, 1] * 1e7 + plist[, 2]
  pid_el <- uid[el.mc[, 1]] * 1e7 + uid[el.mc[, 2]]
  matches <- match(pid_el, pid_plist)

  # initialize so that durat_wks starts with variable distribution
  # currently ptype-agnostic
  if (at == 2) {
    md <- netstats$netmain$durat_wks
    cd <- netstats$netcasl$durat_wks
    me <- netstats$netmain$edges
    ce <- netstats$netcasl$edges
    mw <- me / sum(me, ce)
    cw <- ce / sum(me, ce)

    init_pstart <- weighted.mean(c(md, cd), w = c(mw, cw))
    plist[, "start"] <- rnorm(length(plist[, "start"]), -init_pstart, 100)
  }

  durations <- (at - plist[, "start"])[matches]

  # HIV concordance
  hiv.concord <- rep(0, nrow(el.mc))

  ## record indices for each edge HIV combo
  both.hiv.neg <-
    which(diag.status[el.mc[, 1]] == 0 & diag.status[el.mc[, 2]] == 0)

  both.hiv.pos <-
    which(diag.status[el.mc[, 1]] == 1 & diag.status[el.mc[, 2]] == 1)

  serodisc.hiv <-
    which(diag.status[el.mc[, 1]] != diag.status[el.mc[, 2]])

  ## assign HIV combo to edge
  hiv.concord[both.hiv.neg] <- 1
  hiv.concord[both.hiv.pos] <- 2
  hiv.concord[serodisc.hiv] <- 3

  # Model predictions
  pred_df <- data.table(
    ptype = el.mc[, "ptype"],
    durat_wks = durations,
    race.combo = race.combo,
    age.i = age.i,
    age.j = age.j,
    abs_sqrt_agediff = abs_sqrt_agediff,
    hiv.concord = hiv.concord
  )


  # Predict anal act rates
  # NOTE: exp() used because ai.acts.mod outcome is log(act.rate * 52)
  # ht for elegant simulation from model:
  # https://www.barelysignificant.com/post/glm/

  log.ai.acts <- predict(ai.acts.mod, newdata = pred_df, type = "response")
  log.ai.eps <- rnorm(length(log.ai.acts), 0, sigma(ai.acts.mod))
  ai.rates <- exp(log.ai.acts + log.ai.eps) / 52

  ai.rates <- ai.rates * ai.acts.scale
  ai <- rpois(length(ai.rates), ai.rates)

  # Predict oral act rates
  log.oi.acts <- predict(oi.acts.mod, newdata = pred_df, type = "response")
  log.oi.eps <- rnorm(length(log.oi.acts), 0, sigma(oi.acts.mod))
  oi.rates <- exp(log.oi.acts + log.oi.eps) / 52

  oi.rates <- oi.rates * oi.acts.scale
  oi <- rpois(length(oi.rates), oi.rates)
  el.mc <- cbind(el.mc, durations, ai, oi)

  # Add one-time partnerships
  el.oo <- el[el[, "ptype"] == 3, ]
  ai <- oi <- durations <- rep(1, nrow(el.oo))
  el.oo <- cbind(el.oo, durations, ai, oi)

  # Bind el back together
  el <- rbind(el.mc, el.oo)

  # For AIDS cases with VL above acts.aids.vl, reduce their their anal acts to 0
  p1HIV <- which(el[, "st1"] == 1)
  p1AIDS <-
    stage[el[p1HIV, "p1"]] == 4 & vl[el[p1HIV, "p1"]] >= acts.aids.vl
  el[p1HIV[p1AIDS == TRUE], "ai"] <- 0

  p2HIV <- which(el[, "st2"] == 1)
  p2AIDS <-
    stage[el[p2HIV, "p2"]] == 4 & vl[el[p2HIV, "p2"]] >= acts.aids.vl
  el[p2HIV[p2AIDS == TRUE], "ai"] <- 0

  # Flip order of discordant edges
  disc <- abs(el[, "st1"] - el[, "st2"]) == 1
  disc.st2pos <- which(disc == TRUE & el[, "st2"] == 1)
  el[disc.st2pos, 1:4] <- el[disc.st2pos, c(2, 1, 4, 3)]

  # Remove inactive edges from el (no anal or oral acts)
  el <- el[-which(el[, "ai"] == 0 & el[, "oi"] == 0), ]

  # Save out
  dat$temp$el <- el

  return(dat)
}
