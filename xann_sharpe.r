# Demonstrate ann_sharpe() for calls, puts, and straddles.
#
# The examples compare Black-Scholes option prices at an implied volatility with
# payoff moments under a specified realized volatility. If mu is omitted,
# ann_sharpe() uses the risk-neutral stock drift r - q. The output includes
# option price, payoff mean/sd, profit mean/sd, holding-period Sharpe, and
# annualized Sharpe.

source("option_stats.r")

options(width = 10000)

format_ann_sharpe_demo <- function(x) {
  two_dec <- c("price", "mean_payoff", "sd_payoff", "mean_profit", "sd_profit")
  four_dec <- c("implied_sigma", "mu", "realized_sigma", "T")
  three_dec <- c("sharpe", "ann_sharpe")

  for (col in intersect(two_dec, names(x))) {
    x[[col]] <- sprintf("%.2f", x[[col]])
  }
  for (col in intersect(four_dec, names(x))) {
    x[[col]] <- sprintf("%.4f", x[[col]])
  }
  for (col in intersect(three_dec, names(x))) {
    x[[col]] <- sprintf("%.3f", x[[col]])
  }
  x
}

S0 <- 100
r <- 0.03
q <- 0.00
T <- 1.0
implied_sigma <- 0.20
realized_sigma = 0.25
cat("At-the-money straddle with realized vol =", realized_sigma, "implied vol =", implied_sigma, "\n")
atm_vol_comparison <- do.call(rbind, lapply(c(0.25, 1.0), function(option_T) {
  do.call(rbind, lapply(c("straddle", "call", "put"), function(option_type) {
    ann_sharpe(S0 = S0, K = 100, r = r, q = q, realized_sigma = realized_sigma,
      T = option_T, type = option_type, implied_sigma = implied_sigma)
  }))
}))
print(format_ann_sharpe_demo(atm_vol_comparison), row.names = FALSE)
