
#' @title Position Module
#'
#' @description Module function for establishing sexual role or position in each
#'              act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The sexual role within each act is determined by each nodes "role identity"
#' as exclusively receptive, exclusively insertive, or versatile. This function
#' determines whether the infected or the susceptible partner is the insertive
#' partner for that act. For the first two role identity types, that is
#' deterministic based on identity. For versatile-versatile pairs, this is
#' determined stochastically for each act.
#'
#' @return
#' This function returns the updated discordant edgelist with a \code{ins}
#' attribute for values of whether the infected node is insertive or the
#' susceptible node is insertive for that act.
#'
#' @keywords module msm
#'
#' @export
#'
position_msm <- function(dat, at) {

  ## Variables
  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  # Attributes
  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot

  # Parameters

  ## Process
  p1.role.class <- role.class[al[, 1]]
  p2.role.class <- role.class[al[, 2]]

  ins <- rep(NA, length(p1.role.class))
  ins[which(p1.role.class == "I")] <- 1
  ins[which(p1.role.class == "R")] <- 0
  ins[which(p2.role.class == "I")] <- 0
  ins[which(p2.role.class == "R")] <- 1

  # Versatile MSM
  vv <- which(p1.role.class == "V" & p2.role.class == "V")
  p1.ins.prob <- ins.quot[al[, 1][vv]] /
                 (ins.quot[al[, 1][vv]] + ins.quot[al[, 2][vv]])
  p1.ins <- rbinom(length(vv), 1, p1.ins.prob)
  ins[vv[p1.ins == 1]] <- 1
  ins[vv[p1.ins == 0]] <- 0

  ## Output
  dat$temp$al <- cbind(al, ins)

  return(dat)
}


#' @title Update Role Class in Main and Casual Partnerships
#'
#' @description Module function for updating act class in main and casual
#'              partnerships based on probabilities of transition.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{role.class} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#'
#' @export
#'
update_roleclass_msm <- function(dat, at) {

  role.trans.matrix <- dat$param$role.trans.matrix
  if (sum(colSums(role.trans.matrix) != 1) > 0) {
    stop("Column sums in argument role.trans.matrix must all equal 1.")
  }
  old.role.class <- dat$attr$role.class

  new.role.class <- sapply(1:length(old.role.class),
                           function(x) {
                             sample(c("I", "V", "R"), size = 1,
                                    prob = role.trans.matrix[, 1] *
                                      (old.role.class[x] == "I") +
                                      role.trans.matrix[, 2] *
                                      (old.role.class[x] == "V") +
                                      role.trans.matrix[, 3] *
                                      (old.role.class[x] == "R"))
                           })

  dat$attr$role.class <- new.role.class

  return(dat)
}
