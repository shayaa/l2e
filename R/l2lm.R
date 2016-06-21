#' Fit the L2E Parameters of a Multivariate Linear Model with a Univariate Response
#'
#' @param X a predictor matrix.
#' @param y a response vector.
#' @param init_x an initial estimate of \code{a}, \code{b}, and \code{sig}. Defaults respectively to the mean response \code{mean(y)}, a vector zeros of length \code{ncol(X)}, and the sample standard deviation \code{sqrt(var(y))}
#' @param partial are the residuals modeled by a full density \eqn{\epsilon ~ \mathcal N(0, \sigma^2)} or a partial density \eqn{\epsilon ~ w \mathcal N(0, \sigma^2)} .
#' @param init_w an initial estimate of \code{w}. If \code{partial = FALSE} then w is fixed at 1.
#' @param print_level value of the argument \code{print.level} from the \code{\link[stats]{nlm}} function
#'
#' @return The function \code{l2lm} returns a list of L2E parameters
#' \describe{
#'   \item{a}{the intercept estimate}
#'   \item{b}{a vector of estimated coefficients}
#'   \item{sig}{an estimated standard deviation of the residuals}
#'   \item{w}{the weight of the partial L2E fit; it can be greater than 1}
#'   \item{res}{the residuals of the model}
#' }
#'
#' @examples
#' ## simulate linear model with 20% contamination
#' X = c(3*(1:50) + rnorm(50, 0,20), rnorm(10, 0, 20))
#' y = 1:60
#' l2lm_out = l2lm(X = X, y = y, partial = F)
#' plot(X,y)
#' abline(l2lm_out$a, l2lm_out$b)

l2lm <- function(X, y, init_x, partial  = F,
                 init_w = 1, print_level = 0){

  if (!is.matrix(X)) {
    X <- cbind(X)
    }

  p <- ncol(X)

  if (length(y) != nrow(X)) {
    stop("X and y differ on n")
    }

  l2lm_criteria <- function(x) {
    a <- x[1]
    b <- x[2:(p + 1)]
    sig <- exp(x[p + 2])

    if (partial) {
      w <- exp(x[p + 3])
      } else {
      w <- init_w
      }
    ei <- y - a - X %*% b

    return(w ^ 2 / (2 * sqrt(pi) * sig) -
             2 * w * mean(dnorm(ei, 0, sig)))
  }

  if (missing(init_x)) {
    init_x <- c(mean(y),
             rep(0, p),
             sqrt(var(y)))
    }

  if (length(init_x) != p + 2) {
    stop("init_x error")
  }

  x0 <- c(init_x[1:(p + 1)], log(init_x[p + 2]))

  if (partial) {
    x0 <- c(x0, log(init_w))
    }

  ans <- nlm(l2lm_criteria,
             x0,
             iterlim = 100,
             print.level = print_level)

  a <- ans$est[1]
  b <- ans$est[2:(p + 1)]
  sig <- exp(ans$est[p + 2])
  if (partial) {
    w <- exp(ans$est[p + 3])
    }
  else {
    w <- init_w
    }

  return(list(a = a,
              b = b,
              sig = sig,
              w = w,
              res = c(y - a - X %*% b)))
}
