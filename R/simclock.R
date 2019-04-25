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
#' @references Yang and Rannala (2007) \emph{Inferring speciation times under an
#'   episodic molecular clock.} Systematic Biology, 56:453-466.
#'
#' @return An object of class phylo with branch lengths in substitutions per
#'   site.
#'
#' @author Mario dos Reis
#'
#' @export
relaxed.tree <- function(tree, model, r, s2) {
  tt <- tree
  nb <- length(tt$edge.length)

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
  nb <- tree$Nnode
  nt <- length(tree$tip.label)
  tree$edge.length <- tree$edge.length / 2
  rv <- numeric(nb)
  Sig <- matrix(0, ncol=2, nrow=2)

  for (node in (nt+1):(nt+nb)) {
    dad <- which(tree$edge[,2] == node)
    if (length(dad) == 0) {  # I'm the root!
      ta <- 0  # ancestral time
      ra <- r  # ancestral rate
    }
    else {
      ta <- tree$edge.length[dad]
      ra <- rv[node - nb]
    }

    desc <- which(tree$edge[,1] == node)
    #left <- desc[1]; right <- desc[2]
    tl <- tree$edge.length[desc[1]]
    tr <- tree$edge.length[desc[2]]

    mu <- c(ra - (ta + tl) * s2/2, ra - (ta + tr) * s2/2)
    diag(Sig) <- c(ta + tl, ta + tr) * s2
    Sig[1,2] <- Sig[2,1] <- ta * s2

    rr <- MASS::mvrnorm(1, mu, Sig)
    rv[desc[1]] <- rr[1]; rv[desc[2]] <- rr[2]
  }
  return (exp(rv))
}
