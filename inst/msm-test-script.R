
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV-p")

# Main Test Script ----------------------------------------------------

scr.dir <- "~/Dropbox/Projects/injectable-prep/"
load(file.path(scr.dir, "est/nwstats.rda"))
load(file.path(scr.dir, "est/fit.rda"))
epistats <- readRDS("~/Dropbox/Projects/NetParams/data/artnet.EpiStats.Atlanta.rda")

param <- param_msm(nwstats = st,
                   acts.model = epistats$act.rates,
                   riskh.start = 2,
                   prep.start = 50,
                   prep.start.prob = 0.2,
                   prep.adhr.dist = c(0.089, 0.127, 0.785),
                   prep.adhr.hr = c(0.69, 0.19, 0.02),
                   prep.discont.rate = 1 - (2^(-1/781)),
                   prep.tst.int = 90,
                   prep.risk.int = 182,
                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1)
init <- init_msm(nwstats = st,
                 prev.B = 0.260,
                 prev.W = 0.260)
control <- control_msm(simno = 1,
                       nsteps = 200,
                       nsims = 1,
                       ncores = 1)

sim <- netsim(est, param, init, control)

# df <- as.data.frame(sim)
# names(df)


# Testing/Timing ------------------------------------------------------

# dtemp <- dat
dat <- dtemp
# dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:104) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- test_msm(dat, at)
  dat <- tx_msm(dat, at)
  dat <- progress_msm(dat, at)
  dat <- vl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- trans_msm(dat, at)
  dat <- sti_trans(dat, at)
  dat <- sti_recov(dat, at)
  dat <- sti_tx(dat, at)
  dat <- prevalence_msm(dat, at)
  cat(at, ".", sep = "")
}

length(unique(sapply(dat$attr, length)))
names(dat$attr)
table(dat$attr$prepStat, dat$attr$prepElig)
