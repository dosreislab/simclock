% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simclock.R
\name{relaxed.tree}
\alias{relaxed.tree}
\title{Simulate branch lengths on a phylogeny under a relaxed clock}
\usage{
relaxed.tree(tree, model, r, s2)
}
\arguments{
\item{tree}{an object of class phylo representing a bifurcating phylogeny}

\item{model}{character, the relaxed clock model}

\item{r}{numeric, the mean rate in substitutions per site}

\item{s2}{numeric, the rate "diffusion" parameter for the relaxed clocks}
}
\value{
An object of class phylo with branch lengths in substitutions per
  site.
}
\description{
Simulate branch lengths on a phylogeny under a relaxed clock
}
\details{
The \code{tree} is assummed to be a timetree. Thus, if all your tip
  species are extant, then \code{tree} must be ultrametric. If \code{tree} is
  not ultrametric then it is assummed you have extinct tips. The \code{tree}
  must be rooted and strictly bifurcating.

  The options for \code{model} are "clk", "iln", and "gbm_RY07", for the
  strict clock, the independent log-normal rates, and the geometric-Brownian
  motion rates, respectively (see Rannala and Yang, 2007). If \code{model ==
  "clk"} the branch lengths of \code{tree} are multiplied by \code{r}. If
  \code{model == "iln"} or \code{model == "gbm_RN07"}, \eqn{n = 2*s - 2}
  rates (one for each branch in the \code{s} species phylogeny) are sampled
  from the appropriate distribution. The branch lengths in \code{tree} are
  then multiplied by the corresponding rates.
}
\examples{
require(ape)
par(mfrow=c(2,3))

data(pri10s)
# Simulate using autocorrelated log-normal rates on a primate phylogeny:
tt <- relaxed.tree(pri10s, model="gbm", r=.04e-2, s2=.26e-2)

# The relaxed tree (branch lengths are in substitutions per site):
plot(tt, main="Relaxed primate tree (subs per site)")

# The ultrametric timetree of Primates used in the simulation:
plot(pri10s, main="Primates timetree (Ma)")
axisPhylo() # Time unit in 1 Ma.

# Plot the branch lengths for both trees against each other:
plot(pri10s$edge.length, tt$edge.length,
xlab="Branch lengths (Million years)", ylab="Branch lengths (subs per site)")
abline(0, .04e-2) # the slope is the substitution rate, r

data(flu289s)
# Simulate using independent log-normal rates on an influenza H1N1 phylogeny:
tt2 <- relaxed.tree(flu289s, model="iln", r=.15e-2, s2=.45)

# The relaxed tree of influenza:
plot(tt2, show.tip.label=FALSE, main="Relaxed influenza tree (subs per site)")

# The timetree of Influenza (not ultrametric):
plot(flu289s, show.tip.label=FALSE, main="Influenza H1N1 timetree (y)")
axisPhylo(root.time=1907.35, backward=FALSE) # Time unit in 1 y.

# Plot the branch lengths:
plot(flu289s$edge.length, tt2$edge.length,
xlab="Branch lengths (years)", ylab="Branch lengths (subs per site)")
abline(0, .15e-2) # the slope is the substitution rate, r

}
\references{
Drummond et al. (2006) \emph{Relaxed phylogenetics and dating with
  confidence.} PLoS Biology, 4(5): e88.

  Panchaksaram et al. (2024) \emph{Bayesian selection of relaxed-clock
  models: Distinguishing between independent and autocorrelated rates.}
  bioRxiv.

  Yang and Rannala (2007) \emph{Inferring speciation times under an
  episodic molecular clock.} Systematic Biology, 56: 453-466.
}
\author{
Mario dos Reis
}
