
#' @title Depature Module
#'
#' @description Module function for simulting both general and disease-related
#'              departures, including deaths, among population members.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#' Which nodes have died is determined stochastically for general deaths using
#' draws from a binomial distribution, and deterministically for disease-related
#' deaths after nodes have reach a maximum viral load value set in the
#' \code{vl.fatal} parameter. Additionally, deterministic departures from the
#' population occur at the maximum age implemented in the model.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_msm}} module for bookkeeping
#' purposes.
#'
#' @keywords module msm
#' @export
#'
departure_msm <- function(dat, at) {

  ## General departures
  active <- dat$attr$active
  age <- floor(dat$attr$age)
  race <- dat$attr$race

  asmr <- dat$param$netstats$demog$asmr

  idsElig <- which(active == 1)
  idsEligB <- which(active == 1 & race == 0)
  idsEligW <- which(active == 1 & race == 1)

  rates <- rep(NA, length(idsElig))
  rates[idsEligB] <- asmr[age[idsEligB], "vec.asmr.B"]
  rates[idsEligW] <- asmr[age[idsEligW], "vec.asmr.W"]

  idsDep <- idsElig[rbinom(length(rates), 1, rates) == 1]

  ## HIV-related deaths
  idsDepAIDS <- which(dat$attr$stage == 4 &
                   dat$attr$vl >= dat$param$vl.fatal)

  idsDepAll <- unique(c(idsDep, idsDepAIDS))

  if (length(idsDepAll) > 0) {
    dat$attr$active[idsDepAll] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], idsDepAll)
    }
    dat$attr <- deleteAttr(dat$attr, idsDepAll)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in departures mod")
    }
  }

  ## Summary Output
  dat$epi$dep.gen[at] <- length(idsDep)
  dat$epi$dep.AIDS[at] <- length(idsDepAIDS)

  return(dat)
}


#' @export
#' @rdname departure_msm
deaths_het <- function(dat, at) {

  ### 1. Susceptible Deaths ###

  ## Variables
  male <- dat$attr$male
  age <- dat$attr$age
  cd4Count <- dat$attr$cd4Count

  di.cd4.aids <- dat$param$di.cd4.aids
  ds.exit.age <- dat$param$ds.exit.age

  ## Eligible are: active uninf, pre-death infected, unhealthy old
  idsEligSus <- which((is.na(cd4Count) |
                       cd4Count > di.cd4.aids |
                       (cd4Count <= di.cd4.aids & age > ds.exit.age)))
  nEligSus <- length(idsEligSus)

  # Set age-sex specific rates
  ds.rates <- dat$param$ds.rates
  if (nEligSus > 0) {
    rates <- ds.rates$mrate[100*male[idsEligSus] + age[idsEligSus]]
  }


  ## Process
  nDeathsSus <- 0; idsDeathsSus <- NULL
  if (nEligSus > 0) {
    vecDeathsSus <- which(rbinom(nEligSus, 1, rates) == 1)
    nDeathsSus <- length(vecDeathsSus)
  }


  ## Update Attributes
  if (nDeathsSus > 0) {
    idsDeathsSus <- idsEligSus[vecDeathsSus]
    dat$attr$active[idsDeathsSus] <- 0
  }


  ### 2. Infected Deaths ###

  ## Variables
  active <- dat$attr$active
  di.cd4.rate <- dat$param$di.cd4.rate

  ## Process
  nDeathsInf <- 0; idsDeathsInf <- NULL

  cd4Count <- dat$attr$cd4Count
  stopifnot(length(active) == length(cd4Count))

  idsEligInf <- which(active == 1 & cd4Count <= di.cd4.aids)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecDeathsInf <- which(rbinom(nEligInf, 1, di.cd4.rate) == 1)
    if (length(vecDeathsInf) > 0) {
      idsDeathsInf <- idsEligInf[vecDeathsInf]
      nDeathsInf <- length(idsDeathsInf)
    }
  }

  idsDeathsDet <- which(cd4Count <= 0)
  if (length(idsDeathsDet) > 0) {
    idsDeathsInf <- c(idsDeathsInf, idsDeathsDet)
    nDeathsInf <- nDeathsInf + length(idsDeathsDet)
  }


  ### 3. Update Attributes ###
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
  }

  ## 4. Update Population Structure ##
  inactive <- which(dat$attr$active == 0)
  dat$el[[1]] <- tergmLite::delete_vertices(dat$el[[1]], inactive)
  dat$attr <- deleteAttr(dat$attr, inactive)

  if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
    stop("mismatch between el and attr length in death mod")
  }

  ### 5. Summary Statistics ###
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}