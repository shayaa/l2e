l2lm <- function(X, y, xin, w.opt  = F,
                 win = 1, pl = 1){

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

    if (w.opt) {
      w <- exp(x[p + 3])
      } else {
      w <- win
      }
    ei <- y - a - X %*% b

    return(w ^ 2 / (2 * sqrt(pi) * sig) -
             2 * w * mean(dnorm(ei, 0, sig)))
  }

  if (missing(xin)) {
    xin <- c(mean(y),
             rep(0, p),
             sqrt(var(y)))
    }

  if (length(xin) != p + 2) {
    stop("xin error")
  }

  x0 <- c(xin[1:(p + 1)], log(xin[p + 2]))

  if (w.opt) {
    x0 <- c(x0, log(win))
    }

  ans <- nlm(l2lm_criteria, x0, iterlim = 100, print.level = pl)

  a <- ans$est[1]
  b <- ans$est[2:(p + 1)]
  sig <- exp(ans$est[p + 2])
  if (w.opt) {
    w <- exp(ans$est[p + 3])
    }
  else {
    w <- win
    }

  return(list(a = a, b = b, sig = sig,
              w = w, res = c(y - a - X %*% b)))
}
