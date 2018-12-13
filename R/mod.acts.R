
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
  race <- dat$attr$race
  age <- dat$attr$age
  uid <- dat$attr$uid

  plist <- dat$temp$plist

  # Parameters
  mod <- dat$param$acts.model

  # Construct edgelist
  el.mc <- rbind(dat$el[[1]], dat$el[[2]])
  ptype  <- rep(1:2, times = c(nrow(dat$el[[1]]), nrow(dat$el[[2]])))

  # Base AI rates based on Poisson model for main/casual
  # TODO: switch to ARTnet race encoding
  race.combo <- as.numeric(race[el.mc[, 1]] == "W") +
                as.numeric(race[el.mc[, 2]] == "W")
  comb.age <- age[el.mc[, 1]] + age[el.mc[, 2]]

  # Calculate current partnership durations
  pid_plist <- plist[, 1]*1e7 + plist[, 2]
  pid_el <- uid[el.mc[, 1]]*1e7 + uid[el.mc[, 2]]
  matches <- match(pid_el, pid_plist)
  durations <- (at - plist[, "start"])[matches]

  # Model predictions
  x <- data.frame(ptype = ptype,
                  p_duration = durations,
                  race.combo = race.combo,
                  comb.age = comb.age,
                  city = "Atlanta")
  rates <- unname(predict(mod, newdata = x, type = "response"))/52
  ai <- rpois(length(rates), rates)

  # Add one-time partnerships
  el <- rbind(el.mc, dat$el[[3]])
  ptype <- c(ptype, rep(3, nrow(dat$el[[3]])))
  ai <- c(ai, rep(1, sum(ptype == 3)))

  # Add HIV status
  st1 <- status[el[, 1]]
  st2 <- status[el[, 2]]
  disc <- abs(st1 - st2) == 1
  el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1]
  el <- cbind(el, status[el[, 1]], status[el[, 2]], ptype, ai)
  colnames(el) <- c("p1", "p2", "st1", "st2", "ptype", "ai")

  # Remove inactive edges from el
  el <- el[-which(el[, "ai"] == 0), ]

  # Save out
  dat$temp$el <- el

  return(dat)
}
