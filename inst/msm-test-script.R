
suppressMessages(library("EpiModelHIV"))

# @ TODO:
# eventually import netstats, epistats, and netest into this package
data.dir <- "C:/Users/jason/Documents/Github/egcmsm/egcmsm_artnet"

netstats <- readRDS(file.path(data.dir, "netstats/netstats.Rds"))
est <- readRDS(file.path(data.dir, "netest/netest.Rds"))
epistats <- list()

## Parameters for simplified model without HIV

## @NOTE: Where original parameters had unique values for each race/ethnic group,
##        Other value takes the value assigned to White.

param_xgc <- param_msm(

  # external objects
  netstats = netstats,
  epistats = epistats,

  # demography
  a.rate = 0.00055,   # TODO: tweak this "birth rate" if needed

  # gonorrhea parameters
  rgc.tprob = 0.25,   # TODO: update all transmission probs
  ugc.tprob = 0.25,   # TODO: update all transmission probs
  pgc.tprob = 0.25,   # TODO: update all transmission probs
  rgc.ntx.int = 16.8, # TODO: update all duration probs
  ugc.ntx.int = 16.8, # TODO: update all duration probs
  pgc.ntx.int = 16.8, # TODO: update all duration probs

  # STI testing
  gc.sympt.prob.tx = rep(0.7, 4),   # TODO: update all transmission probs
  gc.asympt.prob.tx = rep(0.2, 4),  # TODO: update all transmission probs

  # HIV testing
  hiv.test.rate = c(
    0.01325,
    0.0125,
    0.0124,
    0.0124
  ), # @ORIG, HIV test rate by race
  hiv.test.late.prob = rep(0.25, 4), # @ORIG, proportion of MSM testing only at late-stage (AIDS)

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
  trans.scale = rep(1.0, 4), # ORIGPARAM
  sti.cond.eff = 0.8,        # TODO: condom efficacy for anal sex
  cond.eff = 0.95,           # ORIGPARAM, condom eff anal HIV transmission
  cond.fail = rep(0.25, 4),
  circ.prob = c(
    0.874,
    0.874,
    0.918,
    0.918
  ) # TODO: Other set to White probability (find alternate value)
)

init_xgc <- init_msm(
  prev.ugc = 0.05,  # TODO: update initialization prevelance
  prev.rgc = 0.05,  # TODO: update initialization prevelance
  prev.pgc = 0.05,  # TODO: update initialization prevelance
  prev.uct = 0,
  prev.rct = 0
)

control_xgc <- control_msm(
  simno = 1,
  nsteps = 10,
  nsims = 1,
  ncores = 1,
  initialize.FUN = initialize_msm,
  aging.FUN = NULL,
  departure.FUN = NULL,
  arrival.FUN = NULL,
  hivtest.FUN = NULL,
  hivtx.FUN = NULL,
  hivprogress.FUN = NULL,
  hivvl.FUN = NULL,
  resim_nets.FUN = NULL,
  acts.FUN = NULL,
  condoms.FUN = NULL,
  position.FUN = NULL,
  prep.FUN = NULL,
  hivtrans.FUN = NULL,
  # stitrans.FUN = stitrans_msm,
  stirecov.FUN = NULL,
  stitx.FUN = NULL,
  #  prev.FUN = prevalence_msm
  # aging.FUN = aging_msm,
  # departure.FUN = departure_msm,
  # arrival.FUN = arrival_msm,
  # # hivtest.FUN = hivtest_msm,
  # # hivtx.FUN = hivtx_msm,
  # # hivprogress.FUN = hivprogress_msm,
  # # hivvl.FUN = hivvl_msm,
  # resim_nets.FUN = simnet_msm,
  # acts.FUN = acts_msm,
  # condoms.FUN = condoms_msm,
  # position.FUN = position_msm,
  # prep.FUN = prep_msm,
  # # hivtrans.FUN = hivtrans_msm,
  # stitrans.FUN = stitrans_msm,
  # stirecov.FUN = stirecov_msm,
  # stitx.FUN = stitx_msm,
  # prev.FUN = prevalence_msm,
  # verbose.FUN = verbose.net
)

sim <- netsim(est, param_xgc, init_xgc, control_xgc)
names(sim$attr[[1]])

# Explore clinical history
par(mar = c(3,3,1,1), mgp = c(2,1,0))
m1 <- sim$temp[[1]]$clin.hist[[1]]
m2 <- sim$temp[[1]]$clin.hist[[2]]
m3 <- sim$temp[[1]]$clin.hist[[3]]
a <- sim$attr[[1]]
h <- which(a$status == 1)


m1[h[1:10], 95:104]
aids <- which(a$stage == 4)
id <- h[58]
plot(m1[id, ], type = "o", ylim = c(1, 7))
data.frame(vl = m1[id, ], stage = m2[id, ], tx = m3[id, ])
a$tt.traj[id]
matplot(t(m1[h[1:500], ]), type = "l", lty = 1, ylim = c(1, 7))


df <- as.data.frame(sim)
names(df)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", mean.smooth = FALSE, ylim = c(0, 1))
plot(sim, y = "num")
plot(sim, y = "dep.gen", mean.smooth = TRUE)
plot(sim, y = "dep.AIDS", mean.smooth = FALSE)
plot(sim, y = "prepCurr")
plot(sim, y = "cc.dx", mean.smooth = FALSE)
plot(sim, y = "cc.linked", mean.smooth = FALSE, ylim = c(0.8, 1))
plot(sim, y = "cc.linked1m", mean.smooth = FALSE)
plot(sim, y = "cc.tx", mean.smooth = FALSE)
plot(sim, y = "cc.tx.any1y", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt1", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt2", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt3", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.dur1y", mean.smooth = FALSE)

plot(sim, y = "hstage.acute", mean.smooth = TRUE)
plot(sim, y = "hstage.chronic", mean.smooth = FALSE)
plot(sim, y = "hstage.aids", mean.smooth = FALSE)

plot(sim, y = "ir100.gc", mean.smooth = FALSE)
plot(sim, y = "ir100.ct", mean.smooth = FALSE)
plot(sim, y = "ir100.sti", mean.smooth = FALSE)
plot(sim, y = "prev.gc", mean.smooth = FALSE)
plot(sim, y = "prev.ct", mean.smooth = FALSE)

plot(sim, type = "formation", network = 1, plots.joined = FALSE)
plot(sim, type = "formation", network = 2, plots.joined = FALSE)
plot(sim, type = "formation", network = 3, plots.joined = FALSE)


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
  dat <- stitrans_msm(dat, at)
  dat <- stirecov_msm(dat, at)
  dat <- stitx_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  verbose.net(dat, "progress", at = at)
}

nrow(dat$temp$plist)
table(dat$temp$plist[, "start"])
table(dat$temp$plist[, "stop"])
head(dat$temp$plist)

plist <- as.data.frame(dat$temp$plist)
pmain <- filter(plist, ptype == 2)
table(pmain$start)
hist(pmain$start)
