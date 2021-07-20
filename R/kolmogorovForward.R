#' Discretize function over space-time grid
#'
#' @param f function to evaluate on grid
#' @param x space variable
#' @param t time variable, can be null to discretize a function over just space.
#'
#' @description {Discretizes a function of space-time over a grid.}
#' @details {Simple double for-loop for space-time functions and one for-loop for functions of a single-variable}
#' @return matrix
gridFunction <- function(f, x, t = NULL)
{
  # If both x and t variable are present
  if(!is.null(t))
  {
    N <- length(t)
    M <- length(x)
    m <- matrix(data = 0, nrow = N, ncol = M)
    for(i in 1:N)
    {
      for(j in 1:M)
      {
        m[i, j] <- f(t[i], x[j])
      }
    }
  } else if(is.null(t)){
    n <- length(x)
    m <- matrix(data = 0, nrow = n)
    for(j in 1:n)
    {
      m[j] <- f(x[j])
    }
  }
  return(m)
}

#' Solve a Feynman-Kac PDE under a Ito diffusion
#'
#' @param region vector defining localized space region, an interval
#' @param maturity the time horizon
#' @param dynamics list of named coefficient functions defining the SDE, see details
#' @param problem list of functions defining the conditional expectation, see details
#' @param variational boolean for optional stopping problems
#' @param resolution the time and space resolution, as a vector, i.e. the number of sub-intervals in the discretization scheme
#'
#' @description {Solves the Kolmogorov-forward PDE for a Feynman-Kac problem with arbitrary
#' discounting rate, running cost, terminal cost, and infinitesimal drift and diffusion coefficient functions.}
#' @details {The argument \code{dynamics} must be a named list containing functions of
#' time and space,\code{(t, x)}, where:
#' \itemize{
#' \item \code{drift} is the infinitesimal drift function,
#' \item \code{diffusion} is the infinitesimal volatility function}
#' and the argument \code{problem} must be a named list containing:
#' \itemize{
#' \item \code{discount} the discounting rate, a function of \code{(t, x)},
#' \item \code{running_cost} the running cost, a function of \code{(t, x)},
#' \item \code{terminal_cost}, the terminal cost, a function of \code{x}}
#' }
#' @return a list containing the time-grid, space-grid, and solution matrix grid.
#' @export feynman_kac
feynman_kac <- function(region, maturity, dynamics, problem, variational = FALSE, resolution = c(100, 100))
{
  # resolution, localization, and step-sizes
  N <- resolution[1]
  M <- resolution[2]
  a <- region[1]
  b <- region[2]
  k <- maturity/N
  h <- (b-a)/M
  # Discretized space and time intervals
  x <- seq(a, b, by = h)
  tt <- seq(0, maturity, by = k)
  # Grids, mu, volat.  1:N+1, 1:M+1; on paper 0:N, 0:M
  grids <- lapply(c(dynamics, problem[-3]), function(y) gridFunction(y, x, maturity-tt))
  names(grids) <- c("m", "v", "r", "f")
  m <- grids$m
  v <- grids$v
  r <- grids$r
  f <- grids$f
  g <- gridFunction(problem[[3]], x)
  # Solution matrix: N+1 x M+1
  u <- matrix(0, nrow = N + 1, ncol = M + 1)
  # Terminal condition is IC in numerical implementation from reversing time:
  u[1, ] <- g
  # Forced boundary conditions for numerical scheme: assumed = I.C.
  u[, 1] <- g[1]
  u[, M + 1] <- g[M + 1]
  alpha <- v^2/(2*h^2)-m/(2*h)
  beta <- -r-v^2/(h^2)
  delta <- v^2/(2*h^2)+m/(2*h)
  for(i in 1:N)
  {

    d <- c(alpha[i, 2]*u[i, 1], rep(0, M-3), delta[i, M]*u[i, M+1])
    ld <- -k*alpha[i, 3:M]
    dd <- 1-k*beta[i, 2:M]
    ud <- -k*delta[i, 2:(M-1)]
    target <- u[i, 2:M] + k*(d + f[i, 2:M])
    u[i+1, 2:M] <- tridiag_solve(ld, dd, ud, target)
    if(variational)
    {
      u[i+1, 2:M] <- pmax(u[i+1, 2:M], g[2:M])
    }
  }
  return(list(t = tt, x = x, u = u))
}


