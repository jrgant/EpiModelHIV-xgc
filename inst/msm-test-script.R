pacman::p_load(
  EpiModelHIV,
  data.table,
  magrittr,
  rms,
  stringr
)

pd_path <- "../egcmsm_artnet"
netstats <- readRDS(paste0(pd_path, "/netstats/netstats.Rds"))
est <- readRDS(paste0(pd_path, "/netest/netest.Rds"))
epistats <- readRDS(paste0(pd_path, "/netstats/epistats.Rds"))

# TODO: Where original parameters had unique values for each race/ethnic group,
#       Other value currently takes the value assigned to White.
param_xgc <- param_msm(
  # external objects
  netstats = netstats,
  epistats = epistats,
  # demography
  a.rate = 0.00052,
  arrival.age = 18,
  # TODO: update all transmission probs (priors for all will be [0, 1])
  u2rgc.tprob = 0.25, # urethral-to-rectal transmission probability
  u2pgc.tprob = 0.25, # urethral-to-pharyngeal transmission probability
  r2ugc.tprob = 0.25, # rectal-to-urethral transmission probability
  p2ugc.tprob = 0.25, # pharyngeal-to-urethral transmission probability
  ## NOTE: Following tprobs used only if the kissing/rimming flags are active
  ## in control_msm.
  r2pgc.tprob = 0.25, # rectal-to-pharyngeal transmission probability
  p2rgc.tprob = 0.25, # pharyngeal-to-rectal transmission probability
  p2pgc.tprob = 0.25, # kissing transmission probability
  # TODO: Update all durations
  rgc.ntx.int = 16.8,
  ugc.ntx.int = 16.8,
  pgc.ntx.int = 16.8,
  # STI testing
  gc.sympt.prob.tx = rep(0.7, 4),   # TODO: update all treatment probs
  gc.asympt.prob.tx = rep(0.2, 4),  # TODO: update all treatment probs
  # HIV testing
  hiv.test.rate = c(
    0.01325,
    0.0125,
    0.0124,
    0.0124
  ), # @ORIG, HIV test rate by race
  # @ORIG, prop. of MSM testing only at late-stage (AIDS)
  hiv.test.late.prob = rep(0.25, 4),
  # HIV treatment parameters
  tt.part.supp = rep(0.2, 4), # ORIGPARAM, partial VLS post ART start
  tt.full.supp = rep(0.4, 4), # ORIGPARAM, full VLS w/ post ART start
  tt.dur.supp = rep(0.4, 4),  # ORIGPARAM, durable VLS post ART start
  tx.init.prob = c(
    0.092,
    0.092,
    0.127,
    0.127
  ), # @ORIG
  tx.halt.part.prob = c(
    0.0102,
    0.0102,
    0.0071,
    0.0071
  ), # @ORIG
  tx.halt.full.rr = rep(0.9, 4), # ORIGPARAM
  tx.halt.dur.rr = rep(0.5, 4),  # ORIGPARAM
  tx.reinit.part.prob = c(
    0.00066,
    0.00066,
    0.00291,
    0.00291
  ), # @ORIG
  tx.reinit.full.rr = rep(1.0, 4), # ORIGPARAM
  tx.reinit.dur.rr = rep(1.0, 4),  # ORIGPARAM
  # scaling parameters
  ai.acts.scale = 1,
  oi.acts.scale = 1,
  kiss.rate.main = 3, # NEWPARAM: Kissing prob. during anal/oral sex, main
  kiss.rate.casl = 3, # NEWPARAM: Kissing prob. during anal/oral sex, casual
  kiss.prob.oo = 0.5, # NEWPARAM: Kissing prob. during anal/oral sex, one-time
  rim.rate.main = 2, # NEWPARAM: Weekly rate of analingus in main partnerships
  rim.rate.casl = 2, # NEWPARAM: Weekly rate of analingus in casual partnerships
  rim.prob.oo = 0.2, # NEWPARAM: Prob. of rimming during one-time sexual contact
  trans.scale = rep(1.0, 4), # ORIGPARAM
  sti.cond.eff = 0.8,        # TODO: condom efficacy for anal sex
  cond.eff = 0.95,           # ORIGPARAM, condom eff anal HIV transmission
  cond.fail = rep(0.25, 4),
  sti.cond.fail = rep(0.2, 4),
  circ.prob = c(
    0.874,
    0.874,
    0.918, # TODO: Other set to White probability (find alternate value)
    0.918
  )
)

init_xgc <- init_msm(
  prev.ugc = 0.1,  # TODO: update initialization prevalence
  prev.rgc = 0.1,  # TODO: update initialization prevalence
  prev.pgc = 0.1   # TODO: update initialization prevalence
)

control_xgc <- control_msm(
  # Computing options
  simno = 1,
  nsteps = 20,
  nsims = 20,
  ncores = 4,
  # Epidemic simulation options
  transRoute_Kissing = TRUE,  # FLAG: Toggle kissing transmission
  transRoute_Rimming = TRUE,  # FLAG: Toggle rimming transmission
  tergmLite = TRUE,  # NOTE Must be set to avoid error thrown by saveout.net()
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
  condoms.FUN = condoms_msm,    # NOTE All act lists are finalized here
  position.FUN = position_msm,
  prep.FUN = prep_msm,
  hivtrans.FUN = hivtrans_msm,
  stirecov.FUN = stirecov_msm,  # TODO Add alternate dists
  stitx.FUN = stitx_msm,        # TODO Add alternate screening/testing protocols
  stitrans.FUN = stitrans_msm,
  prev.FUN = prevalence_msm,
  verbose.FUN = verbose.net
 )

sim <- netsim(est, param_xgc, init_xgc, control_xgc)

saveRDS(sim, "output/dummy_run.Rds")

plot_vec <- function(vec, maint = "") {
  plot(
    sim,
    "epi",
    vec,
    sim.lines = TRUE,
    mean.line = FALSE,
    mean.lwd = 2,
    qnts = FALSE,
    legend = TRUE,
    main = maint
  )
}

races <- c("B", "H", "O", "W")

# GC
plot_vec(paste0("incid.", c("gc", "rgc", "ugc", "pgc")))
plot_vec(paste0("incid.", races, ".rgc"), "Rectal GC")
plot_vec(paste0("incid.", races, ".ugc"), "Urethral GC")
plot_vec(paste0("incid.", races, ".pgc"), "Pharyngeal GC")

plot_vec(
  names(sim$epi) %>% .[. %like% "incid\\.[rup]{1}2(p|r|u)gc"],
  "Transmission pathways (incident infections)"
)

# HIV
plot_vec(paste0("i.num.", races))
plot_vec("i.num")
plot_vec("incid")
plot_vec()

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
  ## dat <- condoms_msm(dat, at)
  ## dat <- position_msm(dat, at)
  ## dat <- prep_msm(dat, at)
  ## dat <- hivtrans_msm(dat, at)
  ## dat <- stitrans_msm(dat, at)
  ## dat <- stirecov_msm(dat, at)
  ## dat <- stitx_msm(dat, at)
  ## dat <- prevalence_msm(dat, at)
  ## verbose.net(dat, "progress", at = at)
}

nrow(dat$temp$plist)
table(dat$temp$plist[, "start"])
table(dat$temp$plist[, "stop"])
head(dat$temp$plist)

plist <- as.data.frame(dat$temp$plist)
pmain <- filter(plist, ptype == 2)
table(pmain$start)
hist(pmain$start)
