
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV-p")

# Main Test Script ----------------------------------------------------

scr.dir <- "~/Dropbox/Projects/NetParams/"
netstats <- readRDS(file.path(scr.dir, "data/artnet.NetStats.Atlanta.rda"))
epistats <- readRDS(file.path(scr.dir, "data/artnet.EpiStats.Atlanta.rda"))
est <- readRDS(file.path(scr.dir, "data/artnet.NetEst.Atlanta.rda"))

file.info(file.path(scr.dir, "data/artnet.NetEst.Atlanta.rda"))$mtime >
  file.info(file.path(scr.dir, "data/artnet.NetStats.Atlanta.rda"))$mtime

param <- param_msm(netstats = netstats,
                   hiv.test.int = c(301, 315),
                   a.rate = 0.00055 / 7,
                   riskh.start = 2,
                   prep.start = 30,
                   prep.start.prob = 0.10,
                   tt.part.supp = c(0.20, 0.20),
                   tt.full.supp = c(0.40, 0.40),
                   tt.dur.supp = c(0.40, 0.40),
                   tx.halt.full.rr = 0.8,
                   tx.halt.dur.rr = 0.1,
                   tx.reinit.full.rr = 2.0,
                   tx.reinit.dur.rr = 5.0,
                   hiv.rgc.rr = 2.5,
                   hiv.ugc.rr = 1.5,
                   hiv.rct.rr = 2.5,
                   hiv.uct.rr = 1.5,
                   hiv.dual.rr = 0.0)
init <- init_msm(init.hiv.mod = epistats$hiv.mod)
control <- control_msm(simno = 1,
                       nsteps = 52 * 25,
                       nsims = 1,
                       ncores = 1)

sim <- netsim(est, param, init, control)

df <- as.data.frame(sim)
names(df)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", mean.smooth = FALSE)
plot(sim, y = "num")
plot(sim, y = "dep.gen", mean.smooth = FALSE)
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

plot(sim, y = "hstage.acute", mean.smooth = FALSE)
plot(sim, y = "hstage.chronic", mean.smooth = FALSE)
plot(sim, y = "hstage.aids", mean.smooth = FALSE)

plot(sim, y = "ir100.gc", mean.smooth = FALSE, ylim = c(0, 10))
plot(sim, y = "ir100.ct", mean.smooth = FALSE, ylim = c(0, 10))

# Testing/Timing ------------------------------------------------------

m <- microbenchmark::microbenchmark(hivvl_msm(dat, at))
print(m, unit = "ms")

dat <- initialize_msm(est, param, init, control, s = 1)

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



