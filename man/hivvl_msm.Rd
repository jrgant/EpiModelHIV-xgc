% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.hivvl.R
\name{hivvl_msm}
\alias{hivvl_msm}
\title{Viral Load Module}
\usage{
hivvl_msm(dat, at)
}
\arguments{
\item{dat}{Master data list object of class \code{dat} containing networks,
individual-level attributes, and summary statistics.}

\item{at}{Current time step.}
}
\value{
This function returns the \code{dat} object with updated \code{vl} attribute.
}
\description{
Module function for updating HIV viral load.
}
\details{
HIV viral load varies over time as a function of time since infection and ART
history. In the absence of ART, VL rises during the acute rising stage and
falls during the acute falling stage, until it reaches a set-point value in
chronic stage infection. VL again rises during AIDS stage disease until the
point of death.

For persons who have ever initated treatment (\code{tt.traj} is \code{3} or
\code{4}), VL changes depending on current ART use in that time step.
Current use is associated with a reduction in VL, with the rates of decline
and nadirs dependent on partial or full suppression levels. Current
non-adherence is associated with an equal level of increase to VL. All persons
who have reached AIDS, regardless of how they arrived, have a similar rate of
VL increase.
}
\keyword{module}
\keyword{msm}
