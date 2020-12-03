
#' @title CD4 Progression Module
#'
#' @description Function used in `hivvl_msm()` to handle CD4 count.
#'

expected_cd4 <- function(method, cd4Count1, cd4Count2,
                         male, age, ageInf,
                         at, time.unit = 7) {

  ## Variables
  timeInf <- (age - ageInf) * (365 / time.unit)
  catAI <- cut(ageInf, breaks = c(0, 30, 40, 50, Inf),
               labels = FALSE, right = FALSE)

  ## Model parameters
  base.male <- 23.53 - 0.76
  base.feml <- base.male + 1.11
  bases <- c(base.feml, base.male)
  ind.bases <- bases[male + 1]

  # Yearly slopes
  slope1 <- -1.49 + 0.34
  slope2 <- slope1 - 0.10
  slope3 <- slope1 - 0.34
  slope4 <- slope1 - 0.63
  slopes <- c(slope1, slope2, slope3, slope4)
  ind.slopes <- slopes[catAI] * (time.unit / 365)

  if (method == "timeto") {
    tt1 <- (sqrt(cd4Count1) - ind.bases)/ind.slopes
    if (!missing(cd4Count2)) {
      tt2 <- (sqrt(cd4Count2) - ind.bases)/ind.slopes
      return(tt2 - tt1)
    } else {
      return(tt1)
    }
  } else {
    if (method == "assign") {
      cd4CountSqrt <- ind.bases + (ind.slopes * timeInf)
      cd4CountSqrt <- pmax(1, cd4CountSqrt)
    }
    if (method == "update") {
      cd4CountSqrt <- sqrt(cd4Count1) + ind.slopes
      cd4CountSqrt[cd4CountSqrt < 1] <- 0
    }
    cd4Count <- cd4CountSqrt ^ 2
    return(cd4Count)
  }

}
