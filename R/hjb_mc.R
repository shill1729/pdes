# # The HJB-PDE fd approx is unstable for non-regular functions
# # and the optimal control is simply infinite for many others.
# # But Monte-Carlo behaves better.
# library(pdes)
# bankroll <- 50
# maturity <- 1
# mu <- 0.1
# rate <- 0.05
# volat <- 0.5
# ic <- function(x) fkpde::indicator(log(x/bankroll) > 0.1)
# ic <- function(x) fkpde::indicator(x > 1.15*bankroll)
# # ic <- function(x) x
# n <- 15000
# hjb_mc <- function(a, maturity, bankroll, mu, rate, volat, ic, n)
# {
#   y <- matrix(0, nrow = length(a))
#   for(i in 1:length(a))
#   {
#     xT <- findistr::rgbm(n, maturity, bankroll, rate+(mu-rate)*a[i], a[i]*volat)
#     y[i] <- mean(ic(xT))
#   }
#   return(y)
# }
#
# hjb_mc_value <- function(bankroll, maturity, mu, rate, volat, ic, n, lb = 0, ub = 10, m = 100)
# {
#   a <- seq(lb, ub, length.out = m)
#   y <- hjb_mc(a, maturity, bankroll, mu, rate, volat, ic, n)
#   a_index <- which.max(y)
#   a_star <- a[a_index]
#   f <- matrix(0, nrow = length(bankroll))
#   for(i in 1:length(bankroll))
#   {
#     xT <- findistr::rgbm(n, maturity, bankroll[i], rate+(mu-rate)*a_star, a_star*volat)
#     f[i] <- mean(ic(xT))
#   }
#   return(f)
# }
# a <- seq(0, 10, length.out = 30)
# y <- hjb_mc(a, maturity, bankroll, mu, rate, volat, ic, n)
# a_index <- which.max(y)
# a_star <- a[a_index]
# b <- seq(1, 2*bankroll, length.out = 30)
# yy <- hjb_mc_value(b, maturity, mu, rate, volat, ic, n)
# par(mfrow = c(1, 2))
# plot(a, y, type = "l", main = "u(t, x, a)")
# abline(v = a_star)
# plot(b, yy, type = "l", xlab = "x", main = "u(t, x)")
# # Output
# print(data.frame(control = a_star, value = hjb_mc(a_star, maturity, bankroll, mu, rate, volat, ic, n)))
