# library(pdes)
# region <- c(-10, 10)
# maturity <- 1
# dynamics <- list(drift = function(t, x) 0,
#                  diffusion  = function(t, x) 1
#                  )
# problem <- list(discount = function(t, x) 0,
#                 running_cost = function(t, x) indicator(x > 0),
#                 terminal_cost = function(x) 0
#                 )
# sol <- feynman_kac(region, maturity, dynamics, problem)
# plot(sol$x, sol$u[100+1, ], type = "l")
# lines(sol$x, sol$u[75, ], lty = "dashed")
# lines(sol$x, sol$u[50, ], lty = "dashed")
# lines(sol$x, sol$u[25, ], lty = "dashed")
#
