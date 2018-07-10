
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at >= dat$param$riskh.start) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < dat$param$prep.start) {
    return(dat)
  }

  # Pull Data ---------------------------------------------------------------

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  race <- dat$attr$race
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepIndic <- dat$attr$prepIndic
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters
  prep.coverage <- dat$param$prep.coverage
  prep.adhr.dist <- dat$param$prep.adhr.dist
  prep.discont.rate <- dat$param$prep.discont.rate
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 &
                        status == 0 &
                        prepStat == 0 &
                        lnt == at)

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind | ind4 >= twind)
  prepIndic[idsIndic] <- 1

  idsEligStart <- intersect(idsIndic, idsEligStart)
  prepElig[idsEligStart] <- 1


  ## Stoppage ------------------------------------------------------------------

  # No indications
  idsNoIndic <- which((ind1 < twind | is.na(ind1)) &
                      (ind2 < twind | is.na(ind2)) &
                      (ind3 < twind | is.na(ind3)) &
                      (ind4 < twind | is.na(ind4)))
  prepIndic[idsNoIndic] <- 0

  # Risk reassessment rules
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at &
                             (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random (memoryless) discontinuation
  idsEligStpRand <- which(active == 1 & prepStat == 1)
  vecStpRand <- rbinom(length(idsEligStpRand), 1, prep.discont.rate)
  idsStpRand <- idsEligStpRand[which(vecStpRand == 1)]

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    idsStart <- ssample(idsEligSt, nStart)
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepIndic <- prepIndic
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected persons
#'              for purpose of PrEP targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Attributes
  n <- length(dat$attr$active)
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.uai.mono)) {
    dat$attr$prep.ind.uai.mono <- rep(NA, n)
    dat$attr$prep.ind.uai.nmain <- rep(NA, n)
    dat$attr$prep.ind.ai.sd <- rep(NA, n)
    dat$attr$prep.ind.sti <- rep(NA, n)
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.uai.mono[part.not.tested.6mo] <- at

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$prep.ind.uai.nmain[uai.nmain] <- at

  ## Condition 3a: AI within known serodiscordant partnerships
  # TODO: remork AI in SD partners list
  # el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  el2.cond3 <- el2[which(el2$st1 == 1 &
                         dat$attr$diag.status[el2$p1] == 1 &
                         el2$ptype %in% 1:2), ]

  ai.sd <- el2.cond3$p2
  dat$attr$prep.ind.ai.sd[ai.sd] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

  return(dat)
}


calc_psens_stats <- function(dat, at) {

  ## Variables

  # params
  prep.timing.lnt <- dat$param$prep.timing.lnt
  prep.indics <- dat$param$prep.indics

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  lnt <- dat$attr$last.neg.test

  prepStat <- dat$attr$prepStat

  if (is.null(dat$attr$prepElig.ever)) {
    dat$attr$prepElig.ever <- rep(0, length(active))
    dat$attr$prepElig.first <- rep(NA, length(active))
    dat$attr$prepElig.last <- rep(NA, length(active))
  }
  prepElig.ever <- dat$attr$prepElig.ever
  prepElig.first <- dat$attr$prepElig.first
  prepElig.last <- dat$attr$prepElig.last

  newBirths <- which(dat$attr$arrival.time == at)
  prepElig.ever[newBirths] <- 0


  # Base eligibility
  if (prep.timing.lnt == TRUE) {
    idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & lnt == at)
  } else {
    idsEligStart <- which(active == 1 & status == 0 & prepStat == 0)
  }

  # Indications
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int

  if (prep.indics == 5) {
    idsEligStart <- intersect(which(ind1 >= twind | ind2 >= twind |
                                      ind3 >= twind | ind4 >= twind),
                              idsEligStart)
  } else if (prep.indics == 1) {
    idsEligStart <- intersect(which(ind1 >= twind),
                              idsEligStart)
  } else if (prep.indics == 2) {
    idsEligStart <- intersect(which(ind2 >= twind),
                              idsEligStart)
  } else if (prep.indics == 3) {
    idsEligStart <- intersect(which(ind3 >= twind),
                              idsEligStart)
  } else if (prep.indics == 4) {
    idsEligStart <- intersect(which(ind4 >= twind),
                              idsEligStart)
  }

  prepElig.ever[idsEligStart] <- 1

  firstPrepElig <- which(is.na(prepElig.first[idsEligStart]))
  prepElig.first[idsEligStart[firstPrepElig]] <- at

  prepElig.last[idsEligStart] <- at


  # Attributes
  dat$attr$prepElig.ever <- prepElig.ever
  dat$attr$prepElig.first <- prepElig.first
  dat$attr$prepElig.last <- prepElig.last


  return(dat)
}
