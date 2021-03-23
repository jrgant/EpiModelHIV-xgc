pacman::p_load(
  EpiModelHIV,
  data.table,
  magrittr,
  rms,
  stringr,
  pscl
)

suppressMessages(library(EpiModelHIV))

an_path   <- "../egcmsm_artnet"
netstats  <- readRDS(paste0(an_path, "/netstats/netstats.Rds"))
est       <- readRDS(paste0(an_path, "/netest/netest.Rds"))
epistats  <- readRDS(paste0(an_path, "/netstats/epistats.Rds"))

param_xgc <- param_msm(
  # external objects
  netstats = netstats,
  epistats = epistats,
  # demography
  a.rate = netstats$demog$mortrate.marginal,
  arrival.age = 18,
  u2rgc.tprob = 0.25, # urethral-to-rectal transmission probability
  u2pgc.tprob = 0.25, # urethral-to-pharyngeal transmission probability
  r2ugc.tprob = 0.25, # rectal-to-urethral transmission probability
  p2ugc.tprob = 0.25, # pharyngeal-to-urethral transmission probability
  ## NOTE: Following tprobs used only if the kissing/rimming flags are active
  ## in control_msm.
  r2pgc.tprob = 0, # rectal-to-pharyngeal transmission probability
  p2rgc.tprob = 0, # pharyngeal-to-rectal transmission probability
  p2pgc.tprob = 0, # kissing transmission probability
  ## Untreated GC durations
  rgc.ntx.int = 15,
  ugc.ntx.int = 2,
  pgc.ntx.int = 12,
  ## Treated GC recovery probabilities c(after 1 week, after 2 weeks, after 3 weeks)
  rgc.tx.recov.pr = c(1 - 0.06, 0.5, 1),
  ugc.tx.recov.pr = c(1 - 0.05, 0.5, 1),
  pgc.tx.recov.pr = c(1 - 0.13, 0.5, 1),
  # STI testing
  gc.sympt.seek.test.prob = 0.71,
  ## NOTE:
  ## Changed the treatment probs to be receiving test conditional on someone's
  ## seeking STI testing. Repeat 3 times, one for each anatomic site.
  ## Tune asymptomatic testing probability to achieve the overall
  ## probability of receiving an STI test at a given anatomic site.
  gc.sympt.prob.test = rep(1, 3),
  ## NOTE: Order of probs: rectal, urethral, pharyngeal
  gc.asympt.prob.test = rep(0.2, 3),  # TODO: calibrate treatment probs
  # HIV treatment parameters
  tt.part.supp = rep(0.2, 4), # ORIGPARAM, partial VLS post ART start
  tt.full.supp = rep(0.4, 4), # ORIGPARAM, full VLS w/ post ART start
  tt.dur.supp = rep(0.4, 4),  # ORIGPARAM, durable VLS post ART start
  tx.init.prob = c(0.092, 0.092, 0.127, 0.127), # @ORIG
  tx.halt.part.prob = c(0.0102, 0.0102, 0.0071, 0.0071), # @ORIG
  tx.halt.full.rr = rep(0.9, 4), # ORIGPARAM
  tx.halt.dur.rr = rep(0.5, 4),  # ORIGPARAM
  tx.reinit.part.prob = c(0.00066, 0.00066, 0.00291, 0.00291), # @ORIG
  tx.reinit.full.rr = rep(1.0, 4), # ORIGPARAM
  tx.reinit.dur.rr = rep(1.0, 4),  # ORIGPARAM
  # scaling parameters
  ai.acts.scale = 1,
  oi.acts.scale = 1,
  kiss.rate.main = 3, # NEWPARAM: Kissing prob. during anal/oral sex, main
  kiss.rate.casl = 0, # NEWPARAM: Kissing prob. during anal/oral sex, casual
  kiss.prob.oo = 0, # NEWPARAM: Kissing prob. during anal/oral sex, one-time
  rim.rate.main = 0, # NEWPARAM: Weekly rate of analingus in main partnerships
  rim.rate.casl = 0, # NEWPARAM: Weekly rate of analingus in casual partnerships
  rim.prob.oo = 0, # NEWPARAM: Prob. of rimming during one-time sexual contact
  trans.scale = rep(1.0, 4), # ORIGPARAM (HIV transmission)
  cdc.sti.int = 12, # NEWPARAM: Regular CDC GC screening interval
  cdc.sti.hr.int = 6, # NEWPARAM: High-risk CDC GC screening interval
  sti.cond.eff = 0.8, # TODO: condom efficacy for anal sex (handle with prior)
  cond.eff = 0.95, # ORIGPARAM, condom eff anal HIV transmission
  cond.fail = rep(0.25, 4),
  # NOTE: Change to 0 to turn off (reflect uncertainty in cond. effect)
  sti.cond.fail = rep(0.2, 4)
)

init_xgc <- init_msm(
  prev.rgc = 0.1,
  prev.ugc = 0.1,
  prev.pgc = 0.1
)

control_xgc <- control_msm(
  # Computing options
  simno = 1,
  nsteps = 50,
  nsims = 10,
  ncores = 4,
  # Epidemic simulation Modules
  initialize.FUN = initialize_msm,
  aging.FUN = aging_msm,
  departure.FUN = departure_msm,
  arrival.FUN = arrival_msm,
  hivtest.FUN = hivtest_msm,
  hivtx.FUN = hivtx_msm,
  hivprogress.FUN = hivprogress_msm,
  hivvl.FUN = hivvl_msm,
  resim_nets.FUN = simnet_msm,
  acts.FUN = acts_msm,
  condoms.FUN = condoms_msm, # NOTE All act lists are finalized here
  position.FUN = position_msm,
  prep.FUN = prep_msm,
  hivtrans.FUN = hivtrans_msm,
  stitrans.FUN = stitrans_msm_rand,
  stirecov.FUN = stirecov_msm,
  stitx.FUN = stitx_msm,
  prev.FUN = prevalence_msm,
  verbose.FUN = verbose.net,
  # Epidemic simulation options
  # NOTE Determines if kissing considered exposure for CDC testing guidelines
  cdcExposureSite_Kissing = FALSE,
  stiScreeningProtocol = "base",
  gcUntreatedRecovDist = "geom",
  tergmLite = TRUE,  # NOTE Must be set to avoid error thrown by saveout.net()
  debug_stitx = FALSE
 )

# Limit number of lines browser prints
options(deparse.max.lines = 5)

sim <- netsim(est, param_xgc, init_xgc, control_xgc)
# saveRDS(sim, "output/dummy_run.Rds")


# Testing/Timing ------------------------------------------------------

m <- microbenchmark::microbenchmark(hivvl_msm(dat, at))
print(m, unit = "ms")

dat <- initialize_msm(est, param_xgc, init_xgc, control_xgc, s = 1)

for (at in 2:200) {
  dat <- aging_msm(dat, at)
  dat <- departure_msm(dat, at)
  dat <- arrival_msm(dat, at)
  dat <- hivtest_msm(dat, at)
  dat <- hivtx_msm(dat, at)
  dat <- hivprogress_msm(dat, at)
  dat <- hivvl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- hivtrans_msm(dat, at)
  dat <- stitrans_msm_rand(dat, at)
  dat <- stirecov_msm(dat, at)
  dat <- stitx_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  verbose.net(dat, "progress", at = at)
}

## nrow(dat$temp$plist)
## table(dat$temp$plist[, "start"])
## table(dat$temp$plist[, "stop"])
## head(dat$temp$plist)

## plist <- as.data.frame(dat$temp$plist)
## pmain <- filter(plist, ptype == 2)
## table(pmain$start)
## hist(pmain$start)
