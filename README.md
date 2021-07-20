
# pdes

<!-- badges: start -->
<!-- badges: end -->

The package pdes provides finite difference solvers for PDEs arising from the Feynman-Kac connection. These PDEs are connected to a stochastic process defined by an SDE which, in general, has some arbitrary drift and diffusion coefficient function of time and a state-variable (e.g. space or price). These same coefficient functions figure in the PDE itself. 

Here is implemented a general solver that is passed user-defined functions that determine the dynamics of the SDE (the two coefficient functions) and the type of conditional expectation problem that figures in the Feynman-Kac connection (a discount, running-cost, and terminal cost function).


## Installation

You can install the latest version from Github with the package devtools:

``` r
devtools::install_github("pdes")
```

## Example

Here we compute the expected time a Brownian motion is positive under one time-unit of observation.

``` r
library(pdes)
region <- c(-10, 10)
maturity <- 1
dynamics <- list(drift = function(t, x) 0,
                 diffusion  = function(t, x) 1
                 )
problem <- list(discount = function(t, x) 0,
                running_cost = function(t, x) indicator(x > 0),
                terminal_cost = function(x) 0
                )
sol <- feynman_kac(region, maturity, dynamics, problem)
plot(sol$x, sol$u[100+1, ], type = "l") # Default resolution is 101x101
lines(sol$x, sol$u[75, ], lty = "dashed")
lines(sol$x, sol$u[50, ], lty = "dashed")
lines(sol$x, sol$u[25, ], lty = "dashed")
```

