# Smoke test for option_stats_nocpp.r.
#
# This script checks that the R-only option statistics/optimizer path can source
# without C++ kernels and run a small long-only plus nonnegative-terminal
# optimization.

source("option_stats_nocpp.r")

S0 <- 100
mu <- 0.03
realized_sigma <- 0.20
T <- 1.0
r <- 0.03
q <- 0.00

K <- c(0, 90, 100, 90, 100)
option_type <- c("call", "call", "call", "put", "put")
implied_sigma <- c(NA, 0.18, 0.18, 0.18, 0.18)

terminal_inputs <- terminal_distribution_inputs(
  terminal_model = "lognormal",
  S0 = S0,
  mu = mu,
  sigma = realized_sigma,
  T = T,
  K = K,
  type = option_type,
  n_scenarios = 101,
  cvar_n_scenarios = 51
)

vol_quotes <- vol_bid_ask_from_quadratic_spread(
  S0 = S0,
  K = K,
  T = T,
  r = r,
  q = q,
  mid_vol = implied_sigma,
  spread_atm = 0.0,
  spread_quad = 0.0
)

option_price <- bs_option_price_vec(
  S0,
  K,
  r,
  q,
  vol_quotes$mid_vol,
  T,
  option_type
)
greeks <- bs_option_greeks_table(
  S0,
  K,
  r,
  q,
  vol_quotes$mid_vol,
  T,
  option_type
)

result <- run_option_optimization_grid(
  m = terminal_inputs$m,
  v = terminal_inputs$v,
  K = K,
  type = option_type,
  option_price = option_price,
  bid_price = option_price,
  ask_price = option_price,
  mid_price = option_price,
  vol_quotes = vol_quotes,
  greeks = greeks,
  S0 = S0,
  r = r,
  q = q,
  T = T,
  expected_payoff = terminal_inputs$expected_payoff,
  cov_payoff = terminal_inputs$cov_payoff,
  dependency_matrix = NULL,
  dependency_diagnostics = NULL,
  payoff_stats = terminal_inputs$payoff_stats,
  payoff_scenarios = terminal_inputs$payoff_scenarios,
  payoff_moment_scenarios = NULL,
  payoff_grid = terminal_inputs$payoff_grid,
  tail_slope = terminal_inputs$tail_slope,
  lower_tail_slope = terminal_inputs$lower_tail_slope,
  cvar_constraints = data.frame(tail_prob = numeric(), max_loss = numeric()),
  cvar_scenarios = terminal_inputs$cvar_scenarios,
  cvar_payoff_scenarios = terminal_inputs$cvar_payoff_scenarios,
  optimization_objectives = c("mean_variance"),
  risk_aversion_utility = c(1.0),
  constraint_modes = c("long_only", "nonnegative_terminal"),
  budget = 1000.0,
  rf_growth = exp(r * T),
  risk_aversion = 4.0,
  min_terminal_wealth = 0.001,
  initial_contracts = rep(0.0, length(K)),
  max_invested_weight = 1.0,
  constrained_optimizer = "nelder_mead",
  constrained_optimizer_max_starts = 2,
  constrained_optimizer_max_iter = 50,
  report_mtm_var = FALSE,
  print_runs = FALSE
)

cat("option_kernels_loaded:", option_kernels_loaded, "\n")
print(
  result$summary_print[, c("constraint_mode", "objective", "obj_value", "optimizer")],
  row.names = FALSE
)
