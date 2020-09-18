#' Simulate branch lengths on a phylogeny under a relaxed clock
#'
#' @param tree an object of class phylo representing a bifurcating phylogeny
#' @param model character, the relaxed clock model
#' @param r numeric, the mean rate in substitutions per site
#' @param s2 numeric, the rate "diffusion" parameter for the relaxed clocks
#'
#' @details The \code{tree} is assummed to be a timetree. Thus, if all your tip
#'   species are extant, then \code{tree} must be ultrametric. If \code{tree} is
#'   not ultrametric then it is assummed you have extinct tips. The \code{tree}
#'   must be rooted and strictly bifurcating.
#'
#'   The options for \code{model} are "clk", "iln", and "gbm_RY07", for the
#'   strict clock, the independent log-normal rates, and the geometric-Brownian
#'   motion rates, respectively (see Rannala and Yang, 2007). If \code{model ==
#'   "clk"} the branch lengths of \code{tree} are multiplied by \code{r}. If
#'   \code{model == "iln"} or \code{model == "gbm_RN07"}, \eqn{n = 2*s - 2}
#'   rates (one for each branch in the \code{s} species phylogeny) are sampled
#'   from the appropriate distribution. The branch lengths in \code{tree} are
#'   then multiplied by the corresponding rates.
#'
#' @references
#'   Drummond et al. (2006) \emph{Relaxed phylogenetics and dating with
#'   confidence.} PLoS Biology, 4(5): e88.
#'
#'   Yang and Rannala (2007) \emph{Inferring speciation times under an
#'   episodic molecular clock.} Systematic Biology, 56: 453-466.
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
#' tt <- relaxed.tree(pri10s, model="gbm", r=.04e-2, s2=.26e-2) # s2=.26
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
#' xlab="Branch lengths (in units of 1 Ma)", ylab="Branch lengths (subs per site)")
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
#' xlab="Branch lengths (in units of 1 y)", ylab="Branch lengths (subs per site)")
#' abline(0, .15e-2) # the slope is the substitution rate, r
#'
#' @author Mario dos Reis
#'
#' @export
# TODO: add indpendent gamma rates model.
relaxed.tree <- function(tree, model, r, s2) {
  tt <- tree
  nb <- length(tt$edge.length)
  model <- match.arg(model, c("clk", "iln", "gbm_RY07"))

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
    rv <- .sim.gmbRY07(tree, r, s2)
    tt$edge.length <- tt$edge.length * rv
    #stop ("gbm_RN07 is not implemented yet")
  }
  else {
    stop (model, " is an unknown clock model")
  }
  return (tt)
}

.sim.gmbRY07 <- function(tree, r, s2) {
  #nb <- tree$Nnode
  nb <- length(tree$edge.length)
  nt <- length(tree$tip.label)
  tree$edge.length <- tree$edge.length / 2
  rv <- numeric(nb)
  Sig <- matrix(0, ncol=2, nrow=2)

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
    #left <- desc[1]; right <- desc[2]
    tl <- tree$edge.length[desc[1]]
    tr <- tree$edge.length[desc[2]]

    mu <- c(ya - (ta + tl) * s2/2, ya - (ta + tr) * s2/2)
    diag(Sig) <- c(ta + tl, ta + tr) * s2
    Sig[1,2] <- Sig[2,1] <- ta * s2

    rr <- MASS::mvrnorm(1, mu, Sig)
    rv[desc[1]] <- rr[1]; rv[desc[2]] <- rr[2]
    #print(c(node, exp(c(ya, rr))))
  }
  return (exp(rv))
}

#' Calculate quantiles of GBM process
#' @export
gbm_RY07q <- function(p, ra, s2, t, log=FALSE) {
  pps <- qnorm(p, mean=log(ra) - t*s2/2, sd=sqrt(s2 * t))
  if (log) {
    return (pps)
  }
  else {
    return (exp(pps))
  }
}
