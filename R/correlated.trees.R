#' Simulate correlated branch lengths among loci on a phylogeny under a relaxed clock
#'
#' @param tree an object of class phylo representing a bifurcating phylogeny
#' @param model character, the relaxed clock model
#' @param r numeric, the mean rate in substitutions per site
#' @param s2 numeric, the rate "diffusion" parameter for the relaxed clocks
#' @param nloci numeric, the number of trees to simulate (one per locus)
#' @param corr, numeric, the correlation of log rates among loci
#'
#' @details A total of \code{nloci} trees are simulated, with the log-rates
#' for branches across trees having correlation \code{corr}.
#'
#' @return A list with two elements: A list of length \code{nloci} of trees of
#'   class phylo with branch lengths in substitutions per site, and a matrix of
#'   branch rates for each locus.
#'
#' @seealso \link{relaxed.tree} to simulate a single tree.
#'
#' @examples
#' require(ape)
#' par(mfrow=c(2,3))
#'
#' data(pri10s)
#' # ILN model:
#' # Simulate using autocorrelated log-normal rates on a primate phylogeny,
#' # with no correlation among three loci:
#' iln0 <- correlated.trees(pri10s, model="iln", r=.04e-2, s2=.1, 3, 0)
#' lapply(iln0$trees, plot)
#' # Repeat with strong correlation among loci:
#' ilnc <- correlated.trees(pri10s, model="iln", r=.04e-2, s2=.1, 3, 0.9)
#' lapply(ilnc$trees, plot)
#'
#' # GBM model:
#' # Simulate using autocorrelated log-normal rates on a primate phylogeny,
#' # with no correlation among three loci:
#' gbm0 <- correlated.trees(pri10s, model="gbm_RY07", r=.04e-2, s2=.26e-2, 3, 0)
#' lapply(gbm0$trees, plot)
#' # Repeat with strong correlation among loci:
#' gbmc <- correlated.trees(pri10s, model="gbm_RY07", r=.04e-2, s2=.26e-2, 3, 0.9)
#' lapply(gbmc$trees, plot)
#'
#'
#'
#' @author Mario dos Reis
#'
#' @export
correlated.trees <- function(tree, model, r, s2, nloci, corr) {
  tt <- tree
  nb <- length(tt$edge.length)
  model <- match.arg(model, c("iln", "gbm_RY07", "gbm0"))

  # construct among loci covariance matrix (p x p):
  R <- matrix(corr, ncol=nloci, nrow=nloci) * s2
  diag(R) <- s2

  # use Eigen decomposition to simulate correlations
  eS <- eigen(R, symmetric = TRUE)
  ev <- eS$values

  # n x p; n: nb; p: nloci
  lrvm <- matrix(nrow=nb, ncol=nloci)

  for (i in 1:nloci) {
    if (model == "iln") {
      lrvm[,i] <- rnorm(nb, 0, 1)
    } else if (model == "gbm_RY07") {
      lrvm[,i] <- .sim.gbmRY07(tree, 1, 1, log=TRUE, drift=TRUE)
      # NOTE: gmbRY07 needs more testing!
    } else if (model == "gbm0") {
      lrvm[,i] <- .sim.gbmRY07(tree, 1, 1, log=TRUE, drift=FALSE)
      # NOTE: gmb0 needs more testing!
    }
  }

  # p x p %*% p x n = p x n
  rvm <- exp( eS$vectors %*% diag(sqrt(pmax(ev, 0)), nloci) %*% t(lrvm) + log(r) - s2/2 )

  tls <- list()
  for (i in 1:nloci) {
    tls[[i]] <- tt
    tls[[i]]$edge.length <- tt$edge.length * rvm[i,]
  }

  return(list(trees=tls, rates=rvm))
}
