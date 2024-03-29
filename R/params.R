
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param netstats Target statistics and related network initialization data from
#'        the standard ARTnet workflow.
#'
#' @param hiv.test.rate HIV testing rates for age groups: [18, 25), [25, 30), [30, 40), [40, 50), [50, 65]
#' @param hiv.test.late.prob Proportion of black/hispanic/other/white MSM who test only
#'        during AIDS stage infection (vector of length 4).
#' @param test.window.int Length of the HIV test window period in weeks.
#' @param tt.part.supp Proportion of black/hispanic/other/white MSM who enter partial viral
#'        suppression category after ART initiation (vector of length 4).
#' @param tt.full.supp Proportion of black/hispanic/other/white MSM who enter full viral
#'        suppression category after ART initiation (vector of length 4).
#' @param tt.dur.supp Proportion of black/hispanic/other/white MSM who enter durable viral
#'        suppression category after ART initiation (vector of length 4).
#'
#' @param tx.init.prob Probability per time step that a black/hispanic/other/white MSM who has
#'        tested positive will initiate treatment (vector of length 4).
#' @param tx.halt.part.prob Probability per time step that black/hispanic/other/white
#'        MSM who have started treatment and assigned to the partial VL suppression
#'        category will stop treatment (vector of length 4).
#' @param tx.halt.full.rr Relative reduction in \code{tx.halt.part.prob} for
#'        black/hispanic/other/white MSM in the full VL suppression category (vector of length 4).
#' @param tx.halt.dur.rr Relative reduction in \code{tx.halt.part.prob} for
#'        black/hispanic/other/white MSM in the durable VL suppression category (vector of length 4).
#' @param tx.reinit.part.prob Probability per time step that a black/hispanic/other/white
#'        MSM who has stopped treatment and assigned to the partial VL suppression
#'        category will restart treatment (vector of length 4).
#' @param tx.reinit.full.rr Relative reduction in \code{tx.reinit.part.prob} for
#'        black/hispanic/other/white MSM in the full VL suppression category (vector of length 4).
#' @param tx.reinit.dur.rr Relative reduction in \code{tx.reinit.part.prob} for
#'        black/hispanic/other/white MSM in the durable VL suppression category (vector of length 4).
#' @param max.time.off.tx.full.int Number of weeks off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.part.int Number of weeks on treatment for a
#'        partial suppressor before onset of AIDS.
#' @param max.time.off.tx.part.int Number of weeks off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of weeks to peak viremia during acute
#'        infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute
#'        infection.
#' @param vl.acute.fall.int Number of weeks from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of weeks to AIDS for a treatment-naive
#'        patient.
#' @param vl.aids.int Duration of AIDS stage infection in weeks.
#' @param vl.aids.peak Maximum viral load during AIDS stage.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param vl.tx.down.slope Number of log10 units that viral load falls per time
#'        step from treatment initiation or re-initiation until the suppression
#'        level is reached (pre-AIDS stages).
#' @param vl.tx.aids.down.slope Number of log10 units that viral load falls per time
#'        step from treatment initiation or re-initiation until the suppression
#'        level is reached (AIDS stage).
#' @param vl.tx.up.slope Number of log10 units that viral load rises per time
#'        step from treatment halting until expected value.
#' @param aids.mr Mortality rate of persons in the AIDS stage who are currently
#'        off ART.
#'
#' @param a.rate Rate at which MSM enter the population.
#' @param arrival.age Age (in years) of new arrivals.
#'
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param trans.scale Relative scalar on base infection probabilities for model
#'        calibration for black/hispanic/other/white men (vector of length 4).
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#'
#' @param cond.eff Relative risk of HIV infection from anal sex when a condom is
#'        used properly (biological efficacy).
#' @param cond.fail Condom failure rates for HIV for black/hispanic/other/white MSM, as a reduction
#'        in the cond.eff parameter (vector of length 4).
#' @param circ.prob Probability that a black/hispanic/other/white new arrival in the population
#'        will be circumcised (vector of length 4).
#' @param act.stopper.prob Probability that an agent with GC symptoms or on Ab treatment for GC stops all sexual activity while either of those conditions are true.
#' @param epistats GLMs for epidemiological parameter from the standard ARTnet workflow.
#' @param acts.aids.vl Viral load level after which sexual act rate goes to zero.
#' @param ai.acts.scale.mc Scalar for main/casual anal sex act rates, for model calibration.
#' @param oi.acts.scale.mc Scalar for main/casual oral sex act rates, for model calibration.
#' @param kiss.rate.main Weekly rate of tongue kissing in main partnerships
#' @param kiss.rate.casl Weekly rate of tongue kissing in casual partnerships
#' @param kiss.prob.oo Probability of tongue kissing during a one-time sexual contact
#' @param rim.rate.main Weekly rate of rimming in main partnerships
#' @param rim.rate.casl Weekly rate of rimming in casual partnerships
#' @param rim.prob.oo Probability of rimming during a one-time sexual contact
#'
#' @param cond.scale Scalar for condom use probability for model calibration.
#'
#' @param riskh.start Time step at which behavioral risk history assessment occurs.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.start.prob Probability of starting PrEP given current indications. Vector of 4 representing Black, Hispanic, Other, White.
#' @param prep.adhr.dist Proportion of men who are low, medium, and high
#'        adherent to PrEP.
#' @param prep.adhr.hr The hazard ratio for infection per act associated with each
#'        level of adherence (from Grant).
#' @param prep.risk.reassess.method Interval for reassessment of risk indications
#'        of active PrEP users, either \code{"none"} for no reassessment,
#'        \code{"inst"} for weekly, or \code{"year"} for year.
#' @param prep.require.lnt If \code{TRUE}, only start on PrEP if current time step is
#'        equal to the last negative test.
#'
#' @param prep.discont.rate Rate of random discontinuation from PrEP. Vector of 4 representing discontinuation rates for Black, Hispanic, Other, and White MSM.
#'
#' @param prep.tst.int Testing interval for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in weeks.
#'
#' @param u2rgc.tprob Probability of gonorrhea transmission per sex act from an infected urethra to an uninfected rectum.
#' @param u2pgc.tprob Probability of gonorrhea transmission per sex act from an infected urethra to an uninfected oropharynx.
#' @param r2ugc.tprob Probability of gonorrhea transmission per sex act from an infected rectum to an uninfected urethra.
#' @param p2ugc.tprob Probability of gonorrhea transmission per sex act from an infected oropharynx to an uninfected urethra.
#' @param r2pgc.tprob Probability of gonorrhea transmission per sex act from an infected rectum to an uninfected oropharynx.
#' @param p2rgc.tprob Probability of gonorrhea transmission per sex act from an infected oropharynx to an uninfected rectum.
#' @param p2pgc.tprob Probability of gonorrhea transmission per sex act from an infected oropharynx to an uninfected oropharynx.
#' @param rgc.tprob Probability of rectal gonorrhea infection per act.
#' @param ugc.tprob Probability of urethral gonorrhea infection per act.
#' @param pgc.tprob Probability of pharyngeal gonorrhea infection per act.

#' @param rgc.sympt.prob Probability of symptoms given infection with rectal
#'        gonorrhea.
#' @param ugc.sympt.prob Probability of symptoms given infection with urethral
#'        gonorrhea.
#' @param pgc.sympt.prob Probability of symptoms given infection with pharyngeal
#'        gonorrhea.
#'
#' @param rgc.ntx.int Median duration in weeks of untreated rectal gonorrhea.
#' @param ugc.ntx.int Median duration in weeks of untreated urethral gonorrhea.
#' @param pgc.ntx.int Median duration in weeks of untreated pharyngeal gonorrhea.
#'
#' @param rgc.tx.recov.pr Rectal GC recovery probability for those on antibiotic treatment in week 1, 2, or 3 of treatment (vector of length 3).
#' @param ugc.tx.recov.pr Urethral GC recovery probability for those on antibiotic treatment in week 1, 2, or 3 of treatment (vector of length 3).
#' @param pgc.tx.recov.pr Pharyngeal GC recovery probability for those on antibiotic treatment in week 1, 2, or 3 of treatment (vector of length 3).
#'
#' @param gc.sympt.prob.test Probability of being tested at symptomatic rectum, urethra, and pharynx.
#' @param gc.asympt.prob.test Probability of being tested at asymptomatic rectum, urethra, and pharynx.
#' @param ugc.sympt.seek.test.prob Probability per week of seeking an STI test due to having STI symptoms.
#' @param rgc.sympt.seek.test.rr Scalar on ugc.sympt.seek.test.prob for those with symptomatic rectal GC.
#' @param pgc.sympt.seek.test.rr Scalar on rgc.sympt.seek.test.prob for those with symptomatic pharyngeal GC.
#' @param cdc.sti.int Regular CDC screening interval (active when stiScreeningProtocol = "cdc" in `control_msm`)
#' @param cdc.sti.hr.int CDC screening interval for "high-risk" MSM (active when stiScreeningProtocol = "cdc" in `control_msm`)
#' @param prep.sti.screen.int Interval in weeks between STI screening at PrEP visits.
#' @param prep.sti.prob.tx Probability of treatment given positive screening during
#'        PrEP visit.
#' @param sti.cond.eff Relative risk of STI infection from anal sex when a condom is
#'        used properly (biological efficacy).
#' @param sti.cond.fail Condom failure rates for STI for black/hispanic/other/white MSM, as
#'        a reduction in the cond.eff parameter (vector of length 4).
#' @param hiv.rgc.rr Relative risk of HIV infection given current rectal gonorrhea.
#' @param hiv.ugc.rr Relative risk of HIV infection given current urethral gonorrhea.
#' @param hiv.dual.rr Additive proportional risk, from 0 to 1, for HIV infection
#'        given dual infection with both gonorrhea and chlamydia.
#'
#' @param ... Additional arguments passed to the function.

#' @return A list object of class \code{param_msm}, which can be passed to EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
#'
param_msm <- function(netstats,

                      # Clinical
                      hiv.test.rate = c(0.0294, 0.0328, 0.0291, 0.0246, 0.0194),
                      hiv.test.late.prob = c(0.160, 0.202, 0.207, 0.222),
                      test.window.int = 21 / 7,
                      tt.part.supp = c(0.20, 0.20, 0.20, 0.20),
                      tt.full.supp = c(0.40, 0.40, 0.40, 0.40),
                      tt.dur.supp = c(0.40, 0.40, 0.40, 0.40),
                      tx.init.prob = c(0.092, 0.092, 0.127, 0.127),
                      tx.halt.part.prob = c(0.0102, 0.0102, 0.0071, 0.0071),
                      tx.halt.full.rr = c(0.9, 0.9, 0.9, 0.9),
                      tx.halt.dur.rr = c(0.5, 0.5, 0.5, 0.9),
                      tx.reinit.part.prob = c(0.00066, 0.00066, 0.00291, 0.00291),
                      tx.reinit.full.rr = c(1.0, 1.0, 1.0, 1.0),
                      tx.reinit.dur.rr = c(1.0, 1.0, 1.0, 1.0),

                      # HIV natural history
                      max.time.off.tx.full.int = 52 * 15,
                      max.time.on.tx.part.int = 52 * 10,
                      max.time.off.tx.part.int = 52 * 10,
                      vl.acute.rise.int = 6.4,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 6.4,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520,
                      vl.aids.int = 104,
                      vl.aids.peak = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      vl.tx.down.slope = 0.25,
                      vl.tx.aids.down.slope = 0.25,
                      vl.tx.up.slope = 0.25,
                      aids.mr = 1/104,

                      # Demographic
                      a.rate = 0.00388, # set to the marginal mortality rate to balance pop. N
                      arrival.age = 18,

                      # HIV transmission prob
                      URAI.prob = 0.0138,
                      UIAI.prob = 0.0011,
                      trans.scale = c(1, 1, 1, 1),
                      acute.rr = 6,
                      circ.rr = 0.4,
                      cond.eff = 0.95,
                      cond.fail = rep(0, 4),
                      circ.prob = c(0.798, 0.428, 0.600, 0.928),

                      # Behavioral (sex acts, condoms)
                      epistats,
                      acts.aids.vl = 5.75,
                      ai.acts.scale.mc = 1,
                      oi.acts.scale.mc = 1,
                      kiss.rate.main = 0,
                      kiss.rate.casl = 0,
                      kiss.prob.oo = 0,
                      rim.rate.main = 0,
                      rim.rate.casl = 0,
                      rim.prob.oo = 0,
                      cond.scale = 1,
                      act.stopper.prob = 0.8,

                      # Per-act GC transmission probabilities
                      u2rgc.tprob = 0.25,
                      u2pgc.tprob = 0.25,
                      r2ugc.tprob = 0.25,
                      p2ugc.tprob = 0.25,
                      r2pgc.tprob = 0,
                      p2rgc.tprob = 0,
                      p2pgc.tprob = 0,

                      # GC intrahost
                      rgc.sympt.prob = 0.16,
                      ugc.sympt.prob = 0.80,
                      pgc.sympt.prob = 0.05,
                      rgc.tx.recov.pr = c(1 - 0.06, 0.5, 1),
                      ugc.tx.recov.pr = c(1 - 0.05, 0.5, 1),
                      pgc.tx.recov.pr = c(1 - 0.13, 0.5, 1),

                      # GC clinic
                      rgc.ntx.int = 16.8,
                      ugc.ntx.int = 16.8,
                      pgc.ntx.int = 16.8,
                      gc.tx.int = 1.4,
                      gc.sympt.prob.test = rep(1, 3),
                      gc.asympt.prob.test = c(0.15, 0.15, 0.15),
                      ugc.sympt.seek.test.prob = 0.72,
                      rgc.sympt.seek.test.rr = 0.5,
                      pgc.sympt.seek.test.rr = 0.5,

                      # Condoms and effect on HIV transmission
                      sti.cond.eff = 0.9,
                      sti.cond.fail = rep(0, 4),
                      hiv.rgc.rr = 2.78,
                      hiv.ugc.rr = 1.73,
                      hiv.dual.rr = 0.2,

                      # PrEP
                      riskh.start = 2,
                      prep.start = 52,
                      prep.start.prob = rep(0.6, 4),
                      prep.adhr.dist = c(0.029, 0.346, 0.625),
                      prep.adhr.hr = c(0.69, 0.19, 0.01),
                      prep.discont.rate = 0.55 / 0.26 * 0.0096 * 0.26 / c(0.262, 0.300, 0.398, 0.424),
                      prep.tst.int = 90 / 7,
                      prep.risk.int = 182 / 7,
                      prep.sti.screen.int = 182 / 7,
                      prep.sti.prob.tx = 1,
                      prep.risk.reassess.method = "year",
                      prep.require.lnt = TRUE,
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param prev.ugc Initial prevalence of urethral gonorrhea.
#' @param prev.rgc Initial prevalence of rectal gonorrhea.
#' @param prev.pgc Initial prevalence of pharyngeal gonorrhea.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_msm}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords msm
#'
#' @export
init_msm <- function(prev.ugc = 0.005,
                     prev.rgc = 0.005,
                     prev.pgc = 0.005,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "init.net"
  return(p)
}


#' @title Epidemic Model Control Settings
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes
#'        if used in conjunction with the \code{EpiModelHPC} package.
#' @param nsims Number of simulations.
#' @param ncores Number of cores per run, if parallelization is used within the
#'        \code{EpiModelHPC} package.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to 1 to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param initialize.FUN Module function to use for initialization of the epidemic
#'        model.
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-realted depatures.
#' @param arrival.FUN Module function for entries into the sexually active population.
#' @param hivtest.FUN Module function for HIV diagnostic disease testing.
#' @param hivtx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param hivprogress.FUN Module function for HIV disease progression.
#' @param hivvl.FUN Module function for HIV viral load evolution.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param hivtrans.FUN Module function to stochastically simulate HIV transmission
#'        over acts given individual and dyadic attributes.
#' @param stitrans.FUN Module function to simulate GC transmission over current
#'        edgelist.
#' @param stirecov.FUN Module function to simulate recovery from GC, heterogeneous
#'        by disease, site, symptoms, and treatment status.
#' @param stitx.FUN Module function to simulate treatment of GC.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param stiScreeningProtocol STI screening/testing protocol. One of "base", "cdc", "symptomatic", "universal".
#' @param cdcExposureSite_Kissing Boolean. Determines whether kissing is considered an eligible sexual exposure when `stiScreeningProtocol = "cdc"`
#' @param gcUntreatedRecovDist One of "geom" or "pois". Determines whether the distribution of untreated GC recovery times is geometric (memoryless) or Poisson (based on current age of infection). Applies to all anatomic sites.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param save.clin.hist Save individual-level clinical history matrices.
#' @param save.transmat Save complete transmission matrix. See help for \code{EpiModel::control.net}.
#' @param truncate.plist Truncate the cumulative partnership list to only include
#'        active partnerships.
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param verbose.int Time step interval for printing model progress. See help for \code{EpiModel::control.net}. Defaults to 1.
#' @param tergmLite Boolean (default = TRUE). Set to avoid error thrown by saveout.net().
#' @param debug_stitx Boolean (default = FALSE). Print interim output from `stitrans_msm_rand()` for debugging purposes.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{control_msm}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
control_msm <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
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
                        condoms.FUN = condoms_msm,
                        position.FUN = position_msm,
                        prep.FUN = prep_msm,
                        hivtrans.FUN = hivtrans_msm,
                        stitx.FUN = stitx_msm,
                        stirecov.FUN = stirecov_msm,
                        stitrans.FUN = stitrans_msm_rand,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose.net,
                        save.nwstats = FALSE,
                        save.clin.hist = FALSE,
                        save.transmat = FALSE,
                        skip.check = FALSE,
                        truncate.plist = TRUE,
                        verbose = TRUE,
                        verbose.int = 1,
                        tergmLite = TRUE,
                        debug_stitx = FALSE,
                        ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  p$save.other <- c("attr", "temp", "el", "p")

  class(p) <- "control.net"
  return(p)
}
