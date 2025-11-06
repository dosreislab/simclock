#' Simulate branch lengths on a phylogeny under a relaxed clock
#'
#' @param tree an object of class phylo representing a bifurcating phylogeny
#' @param model character, the relaxed clock model
#' @param r numeric, the mean rate in substitutions per site
#' @param s2 numeric, the rate "diffusion" parameter for the relaxed clocks
#'
#' @details The \code{tree} is assummed to be a timetree. Thus, if all the tip
#'   species are extant, then \code{tree} must be ultrametric. If \code{tree} is
#'   not ultrametric then it is assummed you have extinct tips. The \code{tree}
#'   must be rooted.
#'
#'   The options for \code{model} are "clk", "iln", "gbm_RY07", and "gbm0", for
#'   the strict clock, the independent log-normal rates, and the
#'   geometric-Brownian motion rates, respectively.
#'
#'   If \code{model == "clk"} the branch lengths of \code{tree} are multiplied
#'   by \code{r}. For the other models, \eqn{n} rates (one for each branch in
#'   the \code{s} species phylogeny) are sampled from the appropriate
#'   distribution. The branch lengths in \code{tree} are then multiplied by the
#'   corresponding rates.
#'
#'   Model "gbm_RY07" is the the one implemented by Rannala and Yang (2007) with
#'   the mean of the rate stabilized across the branches of the tree, i.e., this
#'   version of the GBM process is a martingale. Mean stabilization is achieved
#'   by adding a drift factor to the Brownian process on the log-rate, and thus
#'   the process on the log-rate is no longer a martingale. An alternative
#'   version "gbm0" is provided in which the drift factor is removed. In this
#'   version the process on the log-rate is now a martingale while the process
#'   on the rate is not. The drift factor causes collapsation of the rates to
#'   zero as time goes to infinity.
#'
#' @references
#'   Drummond et al. (2006) \emph{Relaxed phylogenetics and dating with
#'   confidence.} PLoS Biology, 4(5): e88.
#'
#'   Panchaksaram et al. (2024) \emph{Bayesian selection of relaxed-clock
#'   models: Distinguishing between independent and autocorrelated rates.}
#'   Systematic Biology, 74: 323--334.
#'
#'   Rannala and Yang (2007) \emph{Inferring speciation times under an
#'   episodic molecular clock.} Systematic Biology, 56: 453--466.
#'
#' @return An object of class phylo with branch lengths in substitutions per
#'   site.
#'
#' @examples
#' require(ape)
#' par(mfrow=c(2,3))
#'
#' data(pri10s)
#' # Simulate using autocorrelated log-normal rates on a primate phylogeny:
#' tt <- relaxed.tree(pri10s, model="gbm_RY07", r=.04e-2, s2=.26e-2)
#'
#' # The relaxed tree (branch lengths are in substitutions per site):
#' plot(tt, main="Relaxed primate tree (subs per site)")
#'
#' # The ultrametric timetree of Primates used in the simulation:
#' plot(pri10s, main="Primates timetree (Ma)")
#' axisPhylo() # Time unit in 1 Ma.
#'
#' # Plot the branch lengths for both trees against each other:
#' plot(pri10s$edge.length, tt$edge.length,
#' xlab="Branch lengths (Million years)", ylab="Branch lengths (subs per site)")
#' abline(0, .04e-2) # the slope is the substitution rate, r
#'
#' data(flu289s)
#' # Simulate using independent log-normal rates on an influenza H1N1 phylogeny:
#' tt2 <- relaxed.tree(flu289s, model="iln", r=.15e-2, s2=.45)
#'
#' # The relaxed tree of influenza:
#' plot(tt2, show.tip.label=FALSE, main="Relaxed influenza tree (subs per site)")
#'
#' # The timetree of Influenza (not ultrametric):
#' plot(flu289s, show.tip.label=FALSE, main="Influenza H1N1 timetree (y)")
#' axisPhylo(root.time=1907.35, backward=FALSE) # Time unit in 1 y.
#'
#' # Plot the branch lengths:
#' plot(flu289s$edge.length, tt2$edge.length,
#' xlab="Branch lengths (years)", ylab="Branch lengths (subs per site)")
#' abline(0, .15e-2) # the slope is the substitution rate, r
#'
#' @author Mario dos Reis
#'
#' @export
# TODO: add indpendent gamma rates model.
relaxed.tree <- function(tree, model, r, s2, drift) {
  tt <- tree
  nb <- length(tt$edge.length)
  model <- match.arg(model, c("clk", "iln", "gbm_RY07", "gbm0", "gbm_full", "ou"))

  if (!ape::is.rooted(tt)) {
    stop("tree must be rooted")
  }
  if (r < 0) {
    stop ("r must be positive")
  }

  if (model == "clk") {
    tt$edge.length <- tt$edge.length * r
  }
  # r0 = exp(mu + s2/2) -> mu = log(r0) - s2/2
  else if (model == "iln") {
    rv <- rlnorm(nb, meanlog=log(r) - s2/2, sdlog=sqrt(s2))
    tt$edge.length <- tt$edge.length * rv
  }
  else if (model == "gbm_RY07") {
    rv <- .sim.gbmRY07(tree, r, s2, drift=0)
    tt$edge.length <- tt$edge.length * rv
  }
  else if (model == "gbm0") {
    rv <- .sim.gbmRY07(tree, r, s2, drift=0.5*s2)
    tt$edge.length <- tt$edge.length * rv
  }
  else if (model == 'gbm_full') {
    rv <- .sim.gbmRY07(tree, r, s2, drift=drift)
    tt$edge.length <- tt$edge.length * rv
  }
  else if (model == 'ou') {
    #rv <- .sim.ou(tree, r, s2, drift=drift)
    rv <- .sim.ou(tree, r, s2, r_opt, theta)
    tt$edge.length <- tt$edge.length * rv
  }
  return (tt)
}

# Simulate using the GBM rates of Yang and Rannala (2007, Syst. Biol.)
# Guillaume fixed the function to allow simulation on multi-furcating trees.
# Sishuo extended the model to a full GBM with a drift coefficient, representing
# the non-stochastic change in the process. In gbm_full, RY07 and gbm0 are special
# cases with the drift equal to 0 and 0.5*s2, respectively.
.sim.gbmRY07 <- function(tree, r, s2, log=FALSE, drift=0) {

  nb <- length(tree$edge.length)
  nt <- length(tree$tip.label)
  tree$edge.length <- tree$edge.length / 2
  rv <- numeric(nb)

  for (node in (nt+1):(nb+1)) {
    dad <- which(tree$edge[,2] == node)
    if (length(dad) == 0) {  # I'm the root!
      ta <- 0  # ancestral time
      ya <- log(r) # root rate
    }
    else {
      ta <- tree$edge.length[dad]
      ya <- rv[dad]
    }

    desc <- which(tree$edge[,1] == node)
    desc.t <- tree$edge.length[desc]

    # drift == 0 is the YR07 model with stabilised mean (i.e., the gbm
    # process is forced to be a martingale)
    # "r" is added to make it accommodate gbm_full
    # specifically, r = 0, 0.5*s2, user_specified for RY07, gbm0, gbm_full
    mu <- ya + (drift-s2/2) * (ta + desc.t)

    Sig <- matrix(ta * s2, length(desc), length(desc))
    diag(Sig) <- (ta + desc.t) * s2

    rv[desc] <- MASS::mvrnorm(1, mu, Sig)
    #print(c(node, exp(c(ya, rr))))
  }
  if (log) {
    return (rv)
  } else {
    return (exp(rv))
  }
}


#.sim.ou <- function(tree, r, theta, r_opt, s2, log=FALSE) {
.sim.ou <- function(tree, r, s2, r_opt, theta, log=FALSE) {
  # Simulate rate variation under an OU process along a phylogeny
  # tree   : phylo object (ape)
  # r      : root rate
  # theta  : strength of mean reversion
  # r_opt  : exponential of the long-term mean (optimal_rate)
  # s2     : diffusion variance parameter
  # log    : return on log scale if TRUE
  
  nb <- length(tree$edge.length)
  nt <- length(tree$tip.label)
  tree$edge.length <- tree$edge.length / 2
  rv <- numeric(nb)

  for (node in (nt + 1):(nb + 1)) {
    dad <- which(tree$edge[, 2] == node)

    if (length(dad) == 0) {
      ## root case
      tA <- 0
      yA <- log(r)  # rate at the root
    } else {
      tA <- tree$edge.length[dad]
      yA <- rv[dad]
    }

    desc <- which(tree$edge[, 1] == node)
    desc.t <- tree$edge.length[desc]
    n.desc <- length(desc)

    ## Expected means for descendants
    # note that 'mu' is diff from .sim.gbmRY07()
    mu <- log(r_opt)
    means <- (yA - mu) * exp(-theta * (tA + desc.t)) + mu

    ## Cov matrix: Sigma
    Sigma <- matrix(NA, n.desc, n.desc)
    for (i in 1:n.desc) {
      for (j in 1:n.desc) {
        if (i == j) {
          Sigma[i, i] <- (s2 / (2 * theta)) * (1 - exp(-2 * theta * (tA + desc.t[i])))
        } else {
          # Cov btwn descendants
          Sigma[i, j] <- (s2 / (2 * theta)) * 
            (exp(-theta * (desc.t[i] + desc.t[j])) * (1 - exp(-2 * theta * tA)))
        }
      }
    }

    ## descendant trait values
    rv[desc] <- MASS::mvrnorm(1, means, Sigma)
  }

  if (log) return(rv) else return(exp(rv))
}



#' Calculate quantiles of GBM process
#' @export
# TODO: write documentation
gbm_RY07q <- function(p, ra, s2, t, log=FALSE) {
  pps <- qnorm(p, mean=log(ra) - t*s2/2, sd=sqrt(s2 * t))
  if (log) {
    return (pps)
  }
  else {
    return (exp(pps))
  }
}
