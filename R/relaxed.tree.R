#' Simulate branch lengths on a phylogeny under a relaxed clock
#'
#' @param tree An object of class \code{phylo} representing a rooted phylogeny.
#' @param model Character string specifying the relaxed-clock model.
#' @param r Positive numeric value. For the stochastic autocorrelated models,
#'   this is the rate at the root.
#' @param s2 Non-negative numeric value giving the diffusion variance
#'   parameter. For the ILN model, this is the variance of the log-rate.
#' @param log_r_opt Finite numeric value giving the long-term mean log-rate
#'   for the LOU model. The default is \code{log(r)}.
#' @param drift Finite numeric drift coefficient for the full GBM model.
#' @param theta Finite numeric drift target for the GOU model.
#' @param alpha Positive numeric mean-reversion parameter for the LOU and GOU
#'   models.
#'
#' @details
#' The \code{tree} is assumed to be a time tree. Thus, if all tip species are
#' extant, then \code{tree} must be ultrametric. If \code{tree} is not
#' ultrametric, it is assumed that the tree contains extinct or serially
#' sampled tips. The tree must be rooted.
#'
#' The available models are:
#'
#' \itemize{
#'   \item \code{"clk"}: strict molecular clock.
#'   \item \code{"iln"}: independent lognormal branch rates.
#'   \item \code{"gbm_RY07"}: mean-stabilized geometric Brownian motion.
#'   \item \code{"gbm0"}: geometric Brownian motion with zero log-rate drift.
#'   \item \code{"gbm_full"} or \code{"gbm"}: full geometric Brownian motion
#'         with a user-specified drift.
#'   \item \code{"lou"}: Ornstein-Uhlenbeck process on the log-rate.
#'   \item \code{"gou"}: geometric Ornstein-Uhlenbeck process on the rate.
#' }
#'
#' If \code{model == "clk"}, the branch lengths of \code{tree} are multiplied
#' by \code{r}. For the other models, one rate is simulated for each branch,
#' and each temporal branch length is multiplied by its corresponding rate.
#'
#' Model \code{"gbm_RY07"} is the model implemented by Rannala and Yang (2007),
#' with the mean rate stabilized across the branches of the tree. In this
#' version, the rate process is a martingale. Mean stabilization is achieved
#' through a drift correction of \code{-s2/2} on the log-rate.
#'
#' In model \code{"gbm0"}, the drift of the log-rate is zero. Consequently,
#' the log-rate is a martingale, but the rate itself is not.
#'
#' For the LOU model, the log-rate \eqn{X_t = \log R_t} follows
#'
#' \deqn{
#' dX_t = \alpha(\mathrm{log\_r\_opt}-X_t)\,dt
#'        + \sqrt{s^2}\,dB_t.
#' }
#'
#' Thus, \code{log_r_opt} is the stationary mean of the log-rate. The
#' stationary variance of the log-rate is \eqn{s^2/(2\alpha)}.
#'
#' For the GOU model, the rate follows
#'
#' \deqn{
#' dR_t = \alpha(\theta-\log R_t)R_t\,dt
#'        + \sqrt{s^2}R_t\,dB_t.
#' }
#'
#' By Ito's lemma, the log-rate follows
#'
#' \deqn{
#' d\log R_t =
#' \alpha\left(\theta-\frac{s^2}{2\alpha}-\log R_t\right)dt
#' + \sqrt{s^2}\,dB_t.
#' }
#'
#' Therefore, the stationary mean log-rate under the GOU model is
#'
#' \deqn{
#' \mu = \theta-\frac{s^2}{2\alpha}.
#' }
#'
#' If \code{log_r_opt} denotes the stationary mean log-rate, the equivalent
#' GOU parameter is
#'
#' \deqn{
#' \theta=\mathrm{log\_r\_opt}+\frac{s^2}{2\alpha}.
#' }
#'
#' In this parameterization, \code{exp(log_r_opt)} is the stationary median
#' and geometric mean rate, rather than the stationary arithmetic mean rate.
#'
#' @references
#' Drummond et al. (2006) \emph{Relaxed phylogenetics and dating with
#' confidence.} PLoS Biology, 4(5): e88.
#'
#' Panchaksaram et al. (2024) \emph{Bayesian selection of relaxed-clock
#' models: Distinguishing between independent and autocorrelated rates.}
#' Systematic Biology, 74: 323--334.
#'
#' Rannala and Yang (2007) \emph{Inferring speciation times under an
#' episodic molecular clock.} Systematic Biology, 56: 453--466.
#'
#' @return An object of class \code{phylo} with branch lengths in substitutions
#'   per site.
#'
#' @examples
#' \dontrun{
#' require(ape)
#'
#' data(pri10s)
#'
#' # Rannala-Yang geometric Brownian motion:
#' tt <- relaxed.tree(
#'   pri10s,
#'   model = "gbm_RY07",
#'   r = 0.04e-2,
#'   s2 = 0.26e-2
#' )
#'
#' # LOU model:
#' tt_lou <- relaxed.tree(
#'   pri10s,
#'   model = "lou",
#'   r = 0.04e-2,
#'   s2 = 0.26e-2,
#'   log_r_opt = log(0.05e-2),
#'   alpha = 0.5
#' )
#'
#' # Equivalent GOU parameterization:
#' log_r_opt <- log(0.05e-2)
#' alpha <- 0.5
#' s2 <- 0.26e-2
#'
#' theta <- log_r_opt + s2 / (2 * alpha)
#'
#' tt_gou <- relaxed.tree(
#'   pri10s,
#'   model = "gou",
#'   r = 0.04e-2,
#'   s2 = s2,
#'   theta = theta,
#'   alpha = alpha
#' )
#' }
#'
#' @author Mario dos Reis
#'
#' @export
relaxed.tree <- function(
    tree,
    model,
    r,
    s2,
    log_r_opt = log(r),
    drift = 0,
    theta = NULL,
    alpha = NULL
) {
    model <- match.arg(
        model,
        c(
            "clk",
            "iln",
            "gbm_RY07",
            "gbm0",
            "gbm_full",
            "gbm",
            "lou",
            "gou"
        )
    )

    if (!ape::is.rooted(tree)) {
        stop("tree must be rooted", call. = FALSE)
    }

    .validate_numeric_scalar(
        r,
        name = "r",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    if (model %in% c("gbm_full", "gbm")) {
        .validate_numeric_scalar(
            drift,
            name = "drift"
        )
    }

    if (model == "lou") {
        .validate_numeric_scalar(
            log_r_opt,
            name = "log_r_opt"
        )

        .validate_numeric_scalar(
            alpha,
            name = "alpha",
            lower = 0,
            strict = TRUE
        )
    }

    if (model == "gou") {
        .validate_numeric_scalar(
            theta,
            name = "theta"
        )

        .validate_numeric_scalar(
            alpha,
            name = "alpha",
            lower = 0,
            strict = TRUE
        )
    }

    tt <- tree
    nb <- length(tt$edge.length)

    if (model == "clk") {
        tt$edge.length <- tt$edge.length * r
    }

    # If log(R) ~ N(mu, s2), then:
    #
    # E[R] = exp(mu + s2/2).
    #
    # Setting mu = log(r) - s2/2 gives E[R] = r.
    else if (model == "iln") {
        rv <- stats::rlnorm(
            nb,
            meanlog = log(r) - s2 / 2,
            sdlog = sqrt(s2)
        )

        tt$edge.length <- tt$edge.length * rv
    }

    else if (model == "gbm_RY07") {
        rv <- .sim.gbm(
            tree,
            r,
            s2,
            drift = 0
        )

        # Equivalent log-rate drift:
        #
        # log_drift = -s2/2

        tt$edge.length <- tt$edge.length * rv
    }

    else if (model == "gbm0") {
        rv <- .sim.gbm(
            tree,
            r,
            s2,
            drift = 0.5 * s2
        )

        # Equivalent log-rate drift:
        #
        # log_drift = 0

        tt$edge.length <- tt$edge.length * rv
    }

    else if (model %in% c("gbm_full", "gbm")) {
        rv <- .sim.gbm(
            tree,
            r,
            s2,
            drift = drift
        )

        tt$edge.length <- tt$edge.length * rv
    }

    else if (model == "lou") {
        rv <- .sim.lou(
            tree,
            r,
            s2,
            log_r_opt,
            alpha = alpha
        )

        tt$edge.length <- tt$edge.length * rv
    }

    else if (model == "gou") {
        rv <- .sim.gou(
            tree,
            r,
            s2,
            theta,
            alpha = alpha
        )

        tt$edge.length <- tt$edge.length * rv
    }

    return(tt)
}


# -------------------------------------------------------------------------
# Basic validation functions
# -------------------------------------------------------------------------

.validate_numeric_scalar <- function(
    x,
    name,
    lower = NULL,
    strict = FALSE
) {
    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x)) {
        stop(
            sprintf(
                "%s must be a single finite numeric value",
                name
            ),
            call. = FALSE
        )
    }

    if (!is.null(lower)) {
        if (strict && x <= lower) {
            stop(
                sprintf(
                    "%s must be greater than %s",
                    name,
                    lower
                ),
                call. = FALSE
            )
        }

        if (!strict && x < lower) {
            stop(
                sprintf(
                    "%s must be greater than or equal to %s",
                    name,
                    lower
                ),
                call. = FALSE
            )
        }
    }

    invisible(TRUE)
}


.validate_logical_scalar <- function(x, name) {
    if (!is.logical(x) ||
        length(x) != 1L ||
        is.na(x)) {
        stop(
            sprintf(
                "%s must be TRUE or FALSE",
                name
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# -------------------------------------------------------------------------
# Geometric Brownian motion
# -------------------------------------------------------------------------

# Simulate rates using geometric Brownian motion.
#
# The rate process is
#
#   dR_t = drift * R_t dt + sqrt(s2) * R_t dB_t.
#
# Therefore, the log-rate process is
#
#   d log(R_t) = (drift - s2/2) dt + sqrt(s2) dB_t.
#
# Consequently:
#
#   drift = 0
#
# gives the Rannala-Yang mean-stabilized rate process, whereas
#
#   drift = s2/2
#
# gives zero drift on the log-rate.
#
# The log_drift argument is retained as an alternative parameterization.
.sim.gbm <- function(
    tree,
    r,
    s2,
    log = FALSE,
    drift = 0,
    log_drift = NA_real_
) {
    .validate_numeric_scalar(
        r,
        name = "r",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    .validate_numeric_scalar(
        drift,
        name = "drift"
    )

    .validate_logical_scalar(
        log,
        name = "log"
    )

    if (length(log_drift) != 1L) {
        stop(
            "log_drift must be a single numeric value or NA",
            call. = FALSE
        )
    }

    if (!is.na(log_drift)) {
        .validate_numeric_scalar(
            log_drift,
            name = "log_drift"
        )

        drift <- log_drift + 0.5 * s2
    }

    nb <- length(tree$edge.length)
    nt <- length(tree$tip.label)

    # Rates are simulated at branch midpoints.
    tree$edge.length <- tree$edge.length / 2

    rv <- numeric(nb)

    for (node in (nt + 1):(nb + 1)) {
        dad <- which(tree$edge[, 2] == node)

        if (length(dad) == 0) {
            # Root case
            ta <- 0
            ya <- log(r)
        } else {
            ta <- tree$edge.length[dad]
            ya <- rv[dad]
        }

        desc <- which(tree$edge[, 1] == node)
        desc.t <- tree$edge.length[desc]

        # The drift of the log-rate is drift - s2/2.
        mu <- ya + (drift - s2 / 2) * (ta + desc.t)

        Sig <- matrix(
            ta * s2,
            length(desc),
            length(desc)
        )

        diag(Sig) <- (ta + desc.t) * s2

        rv[desc] <- MASS::mvrnorm(
            n = 1,
            mu = mu,
            Sigma = Sig
        )
    }

    if (log) {
        return(rv)
    } else {
        return(exp(rv))
    }
}


# -------------------------------------------------------------------------
# Log-Ornstein-Uhlenbeck model
# -------------------------------------------------------------------------

# Simulate rate variation under an OU process on the log-rate.
#
# Let X_t = log(R_t). The model is
#
#   dX_t = alpha * (log_r_opt - X_t) dt + sqrt(s2) dB_t.
#
# Therefore:
#
#   E[X_infinity]   = log_r_opt
#   Var[X_infinity] = s2/(2*alpha).
#
# exp(log_r_opt) is the stationary median and geometric mean rate.
#
# The stationary arithmetic mean rate is:
#
#   exp(log_r_opt + s2/(4*alpha)).
.sim.lou <- function(
    tree,
    r,
    s2,
    log_r_opt,
    alpha,
    log = FALSE
) {
    .validate_numeric_scalar(
        r,
        name = "r",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    .validate_numeric_scalar(
        log_r_opt,
        name = "log_r_opt"
    )

    .validate_numeric_scalar(
        alpha,
        name = "alpha",
        lower = 0,
        strict = TRUE
    )

    .validate_logical_scalar(
        log,
        name = "log"
    )

    nb <- length(tree$edge.length)
    nt <- length(tree$tip.label)

    # Rates are simulated at branch midpoints.
    tree$edge.length <- tree$edge.length / 2

    rv <- numeric(nb)

    for (node in (nt + 1):(nb + 1)) {
        dad <- which(tree$edge[, 2] == node)

        if (length(dad) == 0) {
            # Root case
            tA <- 0
            yA <- log(r)
        } else {
            tA <- tree$edge.length[dad]
            yA <- rv[dad]
        }

        desc <- which(tree$edge[, 1] == node)
        desc.t <- tree$edge.length[desc]
        n.desc <- length(desc)

        # Stationary mean of the log-rate.
        mu <- log_r_opt

        # Expected log-rates at descendant branch midpoints.
        means <- (yA - mu) *
            exp(-alpha * (tA + desc.t)) +
            mu

        # Covariance matrix for descendant log-rates.
        Sigma <- matrix(
            NA_real_,
            nrow = n.desc,
            ncol = n.desc
        )

        for (i in seq_len(n.desc)) {
            for (j in seq_len(n.desc)) {
                if (i == j) {
                    Sigma[i, i] <-
                        (s2 / (2 * alpha)) *
                        (
                            1 -
                            exp(
                                -2 * alpha *
                                (tA + desc.t[i])
                            )
                        )
                } else {
                    # Covariance between descendants resulting from their
                    # shared evolutionary path.
                    Sigma[i, j] <-
                        (s2 / (2 * alpha)) *
                        exp(
                            -alpha *
                            (desc.t[i] + desc.t[j])
                        ) *
                        (
                            1 -
                            exp(-2 * alpha * tA)
                        )
                }
            }
        }

        rv[desc] <- MASS::mvrnorm(
            n = 1,
            mu = means,
            Sigma = Sigma
        )
    }

    if (log) {
        return(rv)
    } else {
        return(exp(rv))
    }
}


# -------------------------------------------------------------------------
# Geometric Ornstein-Uhlenbeck model
# -------------------------------------------------------------------------

# Simulate rate variation under the geometric OU process
#
#   dR_t = alpha * (theta - log(R_t)) * R_t dt
#          + sqrt(s2) * R_t dB_t.
#
# Let X_t = log(R_t). By Ito's lemma:
#
#   dX_t =
#       alpha * (theta - s2/(2*alpha) - X_t) dt
#       + sqrt(s2) dB_t.
#
# Thus, the stationary mean of the log-rate is
#
#   mu = theta - s2/(2*alpha).
#
# If log_r_opt is the desired stationary mean log-rate, then
#
#   theta = log_r_opt + s2/(2*alpha).
#
# In that case, exp(log_r_opt) is the stationary median and geometric
# mean rate, not the stationary arithmetic mean rate.
#
# The stationary arithmetic mean rate under the GOU parameterization is
#
#   exp(theta - s2/(4*alpha)).
.sim.gou <- function(
    tree,
    r,
    s2,
    theta,
    alpha,
    log = FALSE
) {
    .validate_numeric_scalar(
        r,
        name = "r",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    .validate_numeric_scalar(
        theta,
        name = "theta"
    )

    .validate_numeric_scalar(
        alpha,
        name = "alpha",
        lower = 0,
        strict = TRUE
    )

    .validate_logical_scalar(
        log,
        name = "log"
    )

    nb <- length(tree$edge.length)
    nt <- length(tree$tip.label)

    # Rates are simulated at branch midpoints.
    tree$edge.length <- tree$edge.length / 2

    rv <- numeric(nb)

    for (node in (nt + 1):(nb + 1)) {
        dad <- which(tree$edge[, 2] == node)

        if (length(dad) == 0) {
            # Root case
            tA <- 0
            yA <- log(r)
        } else {
            tA <- tree$edge.length[dad]
            yA <- rv[dad]
        }

        desc <- which(tree$edge[, 1] == node)
        desc.t <- tree$edge.length[desc]
        n.desc <- length(desc)

        # Stationary mean of the log-rate obtained using Ito's lemma.
        mu <- theta - s2 / (2 * alpha)

        # Expected log-rates at descendant branch midpoints.
        means <- (yA - mu) *
            exp(-alpha * (tA + desc.t)) +
            mu

        # Covariance matrix for descendant log-rates.
        Sigma <- matrix(
            NA_real_,
            nrow = n.desc,
            ncol = n.desc
        )

        for (i in seq_len(n.desc)) {
            for (j in seq_len(n.desc)) {
                if (i == j) {
                    Sigma[i, i] <-
                        (s2 / (2 * alpha)) *
                        (
                            1 -
                            exp(
                                -2 * alpha *
                                (tA + desc.t[i])
                            )
                        )
                } else {
                    # Covariance between descendants resulting from their
                    # shared evolutionary path.
                    Sigma[i, j] <-
                        (s2 / (2 * alpha)) *
                        exp(
                            -alpha *
                            (desc.t[i] + desc.t[j])
                        ) *
                        (
                            1 -
                            exp(-2 * alpha * tA)
                        )
                }
            }
        }

        rv[desc] <- MASS::mvrnorm(
            n = 1,
            mu = means,
            Sigma = Sigma
        )
    }

    if (log) {
        return(rv)
    } else {
        return(exp(rv))
    }
}


# -------------------------------------------------------------------------
# Previous GBM function
# -------------------------------------------------------------------------

# Previous GBM interface retained for backward compatibility.
#
# If drift = TRUE, the mean-stabilized Rannala-Yang model is used.
#
# If drift = FALSE, the log-rate has zero drift.
.sim.gbmRY07 <- function(
    tree,
    r,
    s2,
    log = FALSE,
    drift
) {
    .validate_numeric_scalar(
        r,
        name = "r",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    .validate_logical_scalar(
        log,
        name = "log"
    )

    .validate_logical_scalar(
        drift,
        name = "drift"
    )

    nb <- length(tree$edge.length)
    nt <- length(tree$tip.label)

    tree$edge.length <- tree$edge.length / 2

    rv <- numeric(nb)

    for (node in (nt + 1):(nb + 1)) {
        dad <- which(tree$edge[, 2] == node)

        if (length(dad) == 0) {
            # Root case
            ta <- 0
            ya <- log(r)
        } else {
            ta <- tree$edge.length[dad]
            ya <- rv[dad]
        }

        desc <- which(tree$edge[, 1] == node)
        desc.t <- tree$edge.length[desc]

        if (drift) {
            mu <- ya - (ta + desc.t) * s2 / 2
        } else {
            mu <- ya
        }

        Sig <- matrix(
            ta * s2,
            length(desc),
            length(desc)
        )

        diag(Sig) <- (ta + desc.t) * s2

        rv[desc] <- MASS::mvrnorm(
            n = 1,
            mu = mu,
            Sigma = Sig
        )
    }

    if (log) {
        return(rv)
    } else {
        return(exp(rv))
    }
}


# -------------------------------------------------------------------------
# GBM quantiles
# -------------------------------------------------------------------------

#' Calculate quantiles of the Rannala-Yang GBM process
#'
#' @param p Numeric vector of probabilities between zero and one.
#' @param ra Positive ancestral rate.
#' @param s2 Non-negative diffusion variance parameter.
#' @param t Non-negative elapsed time.
#' @param log Logical. If \code{TRUE}, return quantiles on the log scale.
#'
#' @return A numeric vector of quantiles.
#'
#' @export
gbm_RY07q <- function(
    p,
    ra,
    s2,
    t,
    log = FALSE
) {
    if (!is.numeric(p) ||
        length(p) < 1L ||
        anyNA(p) ||
        any(!is.finite(p)) ||
        any(p < 0 | p > 1)) {
        stop(
            "p must contain probabilities between zero and one",
            call. = FALSE
        )
    }

    .validate_numeric_scalar(
        ra,
        name = "ra",
        lower = 0,
        strict = TRUE
    )

    .validate_numeric_scalar(
        s2,
        name = "s2",
        lower = 0,
        strict = FALSE
    )

    .validate_numeric_scalar(
        t,
        name = "t",
        lower = 0,
        strict = FALSE
    )

    .validate_logical_scalar(
        log,
        name = "log"
    )

    pps <- stats::qnorm(
        p,
        mean = log(ra) - t * s2 / 2,
        sd = sqrt(s2 * t)
    )

    if (log) {
        return(pps)
    } else {
        return(exp(pps))
    }
}
