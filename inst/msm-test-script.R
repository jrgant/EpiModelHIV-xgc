
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
                   riskh.start = 500,
                   prep.start = 500,
                   prep.start.prob = 0.2,
                   prep.adhr.dist = c(0.089, 0.127, 0.785),
                   prep.adhr.hr = c(0.69, 0.19, 0.02),
                   prep.discont.rate = 1 - (2^(-1/781)),
                   prep.tst.int = 90,
                   prep.risk.int = 182,
                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,
                   a.rate = 0.0004 / 7)
init <- init_msm(init.hiv.mod = epistats$hiv.mod)
control <- control_msm(simno = 1,
                       nsteps = 500,
                       nsims = 1,
                       ncores = 1)

# sim <- netsim(est, param, init, control)

# df <- as.data.frame(sim)
# names(df)


# Testing/Timing ------------------------------------------------------

library(microbenchmark)
m <- microbenchmark(hivvl_msm(dat, at))
print(m, unit = "ms")

dat <- initialize_msm(est, param, init, control, s = 1) # check

# for (at in 2:104) {
  at = 2
  dat <- aging_msm(dat, at)        # check
  dat <- departure_msm(dat, at)    # check
  dat <- arrival_msm(dat, at)      # check
  dat <- hivtest_msm(dat, at)      # check
  dat <- hivtx_msm(dat, at)        # check
  dat <- hivprogress_msm(dat, at)  # check
  dat <- hivvl_msm(dat, at)        # check
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
  cat(at, ".", sep = "")
# }

nrow(dat$temp$plist)
table(dat$temp$plist[, "start"])
table(dat$temp$plist[, "stop"])
head(dat$temp$plist)

