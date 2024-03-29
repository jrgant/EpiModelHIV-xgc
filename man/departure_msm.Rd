% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.departure.R
\name{departure_msm}
\alias{departure_msm}
\title{Depature Module}
\usage{
departure_msm(dat, at)
}
\arguments{
\item{dat}{Master data list object of class \code{dat} containing networks,
individual-level attributes, and summary statistics.}

\item{at}{Current time step.}
}
\value{
This function returns the updated \code{dat} object accounting for deaths.
The deaths are deactivated from the main and casual networks, as those are in
\code{networkDynamic} class objects; dead nodes are not deleted from the
instant network until the \code{\link{simnet_msm}} module for bookkeeping
purposes.
}
\description{
Module function for simulting both general and disease-related
             departures, including deaths, among population members.
}
\details{
Deaths are divided into two categories: general deaths, for which demographic
data on age-specific mortality rates applies; and disease-related diseases,
for which the rate of death is a function of progression to end-stage AIDS.
}
\keyword{module}
\keyword{msm}
