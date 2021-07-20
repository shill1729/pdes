# # Finite-difference approximation to HJB-PDE for optimal control
# # of a portfolio trading a GBM stock and a risk-free cash account.
# library(pdes)
# bankroll <- 100
# maturity <- 1/252
# mu <- 0.1
# rate <- 0.05
# volat <- 0.5
# ic <- function(x) log(x)
# rc <- function(t, x) 0
# N <- 50
# M <- 25
# epsilon <- 10^-6
# lambda <- abs(mu-rate)/volat
# k <- maturity/N
# h <- 2*(bankroll-epsilon)/M
# # Localization of space to (epsilon, 2*bankroll-epsilon) and discrete partition:
# x <- seq(epsilon, 2*bankroll-epsilon, length.out = M+1)
# # Initializing solution matrix with IC and BC
# u <- matrix(0, nrow = N+1, ncol = M+1)
# # Initial conditions:
# u[1, ] <- ic(x)
# # Boundary conditions: assumed to be = to IC
# u[, 1] <- ic(x[1])
# u[, M+1] <- ic(x[M+1])
# # Coefficients of lower and upper diagonal for the tridiagonal matrix
# alpha <- k*rate*x/(2*h)
# # Auxillary vector: c(alpha_1 u_0, 0,...,0, alpha_M-1 u_M)
# d <- c(alpha[2]*u[1, 1], rep(0, M-2), -alpha[M]*u[1, M+1])
# alphaV <- alpha[3:M]
# betaV <- rep(1, M-1)
# deltaV <- -alpha[2:(M-1)]
# for(i in 1:N)
# {
#   # The recursive RHS vector of the solution at the previous time step
#   bb <- matrix(0, nrow = M-1)
#   for (j in 2:M)
#   {
#     # Non-linear term: square of first derivative divided by second derivative
#     u_x <- u[i, j + 1] - u[i, j - 1]
#     u_xx <- (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1])
#     nlq <- (u_x)^2 / u_xx
#     if(is.nan(nlq))
#     {
#       stop("The FD-approx to the non-linear term is NaN")
#     }
#     bb[j - 1] <- u[i, j] - k*(lambda^2)*nlq/8.0 - d[j - 1]+rc((i-1)*k, x[j])
#   }
#   # Solve the linear equation using any method; Thomas algorithm works best for tridiagonal systems
#   u[i+1, 2:M] <- tridiag_solve(alphaV, betaV, deltaV, bb)
# }
#
# # Exact solution for log-growth to verify convergence
# exactLog <- function(t, x)
# {
#   log(x)+(maturity-t)*(rate+0.5*lambda^2)
# }
#
# # Optimal control
# index <- M/2+1
# value <- u[N+1, index]
# vx <- (u[N+1, index+1]-u[N+1, index-1])/(2*h)
# vxx <- (u[N+1, index+1]-2*u[N+1, index]+u[N+1, index-1])/(h^2)
# control <- -((mu-rate)*vx)/(x[index]*vxx*volat^2)
# exactControl <- (mu-rate)/volat^2
# print(data.frame(control, exactControl))
#
# # Verify known case: log-growth
# # plot(x, exactLog(0, x), type = "l")
# # lines(x, u[N+1, ], col = "blue", lty = "dashed")
#
# # Plot entire solution
# tt <- seq(0, maturity, by = k)
# plot3D::persp3D(x = tt, y = x, z = u)
# # plot(x, u[N+1, ], type = "l")
# # print(u[N+1, index])
#
