# Multi-stock single-tenor option portfolio prototype.
#
# Terminal stock prices are simulated from a multivariate lognormal model with
# a user-specified covariance matrix. Payoff scenarios are then passed to the
# existing optimizer. The default prototype compares long-only and long/short
# constraints so cheap/rich implied vols across stocks can be bought and sold.

source("option_stats.r")
options(width = 10000)

underlyings <- c("AAA", "BBB")
S0 <- c(AAA = 100, BBB = 80)
mu <- c(AAA = 0.06, BBB = 0.04)
realized_vol <- c(AAA = 0.25, BBB = 0.30)
correlation <- matrix(
  c(1.00, 0.45,
    0.45, 1.00),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(underlyings, underlyings)
)
realized_cov <- diag(realized_vol) %*% correlation %*% diag(realized_vol)
dimnames(realized_cov) <- list(underlyings, underlyings)

T <- 1.0
r <- 0.03
q <- c(AAA = 0.00, BBB = 0.01)

option_chains <- list(
  AAA = data.frame(
    type = rep(c("call", "put"), each = 3),
    strike = rep(c(90, 100, 110), 2),
    implied_vol = rep(c(0.22, 0.20, 0.22), 2)
  ),
  BBB = data.frame(
    type = rep(c("call", "put"), each = 3),
    strike = rep(c(70, 80, 90), 2),
    implied_vol = rep(c(0.36, 0.34, 0.36), 2)
  )
)

n_scenarios <- 1001
cvar_n_scenarios <- 1001
scenario_seed <- 123
cvar_seed <- 456
report_mtm_var <- TRUE
var_horizon_days <- 1.0
var_trading_days <- 252.0
var_n_scenarios <- 2001
var_return_model <- "normal_log"
var_terminal_floor <- 0.0
var_conf_levels <- c(0.95, 0.99)
var_seed <- 789

optimization_objectives <- c("mean_variance", "sharpe")
risk_aversion_utility <- c(1.0)
constraint_modes <- c("long_only", "long_short")
risk_aversion <- 4.0
budget <- 1000.0
max_invested_weight <- 1.0
force_full_investment <- TRUE
initial_contracts <- NULL

use_bid_ask_vols <- TRUE
vol_spread_atm <- 0.01
vol_spread_quad <- 0.10
min_bid_vol <- 0.001
stock_bid_ask_spread <- 0.01

print_zero_weight_options <- FALSE
zero_weight_tol <- 1e-10
write_summary_csv <- FALSE
summary_csv_file <- "multistock_objective_summary.csv"

cat("Multi-stock option optimization\n")
cat("underlyings:", underlyings, "\n")
cat("S0:", S0, "\n")
cat("mu:", mu, "\n")
cat("realized_vol:", realized_vol, "\n")
cat("T:", T, "\n")
cat("r:", r, "\n")
cat("q:", q, "\n")
cat("report_mtm_var:", report_mtm_var, "\n")
cat("var_horizon_days:", var_horizon_days, "\n")
cat("var_n_scenarios:", var_n_scenarios, "\n")
cat("var_conf_levels:", var_conf_levels, "\n")
cat("realized_cov:\n")
print(realized_cov)
cat("\n")

stock_rows <- do.call(rbind, lapply(underlyings, function(sym) {
  data.frame(
    underlying = sym,
    type = "call",
    strike = 0.0,
    implied_vol = NA_real_
  )
}))
option_rows <- do.call(rbind, lapply(underlyings, function(sym) {
  chain <- option_chains[[sym]]
  data.frame(
    underlying = sym,
    type = chain$type,
    strike = chain$strike,
    implied_vol = chain$implied_vol
  )
}))
instrument_table <- rbind(stock_rows, option_rows)
instrument_underlying <- instrument_table$underlying
K <- instrument_table$strike
option_type <- instrument_table$type
instrument_implied_sigma <- instrument_table$implied_vol

if (is.null(initial_contracts)) {
  initial_contracts <- rep(0.0, nrow(instrument_table))
}

terminal_prices <- multivariate_lognormal_terminal_prices(
  S0 = S0,
  mu = mu,
  cov_matrix = realized_cov,
  T = T,
  n_scenarios = n_scenarios,
  seed = scenario_seed
)
colnames(terminal_prices) <- underlyings
cvar_prices <- multivariate_lognormal_terminal_prices(
  S0 = S0,
  mu = mu,
  cov_matrix = realized_cov,
  T = T,
  n_scenarios = cvar_n_scenarios,
  seed = cvar_seed
)
colnames(cvar_prices) <- underlyings

horizon_inputs <- multistock_horizon_inputs(
  terminal_prices = terminal_prices,
  cvar_prices = cvar_prices,
  underlying = instrument_underlying,
  K = K,
  type = option_type,
  underlying_names = underlyings
)

vol_quotes <- data.frame(
  mid_vol = instrument_implied_sigma,
  bid_vol = instrument_implied_sigma,
  ask_vol = instrument_implied_sigma,
  vol_spread = rep(0.0, length(K))
)
for (sym in underlyings) {
  idx <- which(instrument_underlying == sym & K > 0)
  one <- vol_bid_ask_from_quadratic_spread(
    S0 = S0[[sym]],
    K = K[idx],
    T = T,
    r = r,
    q = q[[sym]],
    mid_vol = instrument_implied_sigma[idx],
    spread_atm = if (use_bid_ask_vols) vol_spread_atm else 0.0,
    spread_quad = if (use_bid_ask_vols) vol_spread_quad else 0.0,
    min_bid_vol = min_bid_vol
  )
  vol_quotes[idx, ] <- one
}

mid_price <- multistock_price_vec(S0, instrument_underlying, K, option_type, r, q, vol_quotes$mid_vol, T, underlyings)
bid_price <- multistock_price_vec(S0, instrument_underlying, K, option_type, r, q, vol_quotes$bid_vol, T, underlyings)
ask_price <- multistock_price_vec(S0, instrument_underlying, K, option_type, r, q, vol_quotes$ask_vol, T, underlyings)
stock_idx <- which(K == 0 & option_type == "call")
if (length(stock_idx) > 0) {
  for (idx in stock_idx) {
    sym <- instrument_underlying[idx]
    mid_price[idx] <- S0[[sym]]
    bid_price[idx] <- pmax(S0[[sym]] - 0.5 * stock_bid_ask_spread, 0.0)
    ask_price[idx] <- S0[[sym]] + 0.5 * stock_bid_ask_spread
  }
}
option_price <- ask_price

greeks <- multistock_greeks_table(S0, instrument_underlying, K, option_type, r, q, vol_quotes$mid_vol, T, underlyings)

mtm_var_inputs <- if (report_mtm_var) {
  list(
    mode = "multistock",
    current_mid_price = mid_price,
    S0 = S0,
    underlying = instrument_underlying,
    r = r,
    q = q,
    mid_vol = vol_quotes$mid_vol,
    T = T,
    var_mu = mu,
    var_cov_matrix = realized_cov,
    underlying_names = underlyings,
    horizon_days = var_horizon_days,
    trading_days = var_trading_days,
    n_scenarios = var_n_scenarios,
    return_model = var_return_model,
    terminal_floor = var_terminal_floor,
    conf_levels = var_conf_levels,
    seed = var_seed
  )
} else {
  NULL
}

cat("Instrument universe\n")
print(
  data.frame(
    id = seq_along(K),
    underlying = instrument_underlying,
    type = option_type,
    strike = K,
    mid_vol = round(vol_quotes$mid_vol, 4),
    mid_price = round(mid_price, 4),
    bid_price = round(bid_price, 4),
    ask_price = round(ask_price, 4)
  ),
  row.names = FALSE
)
cat("\n")

state_grid <- rbind(
  apply(terminal_prices, 2, quantile, probs = 0.01),
  apply(terminal_prices, 2, quantile, probs = 0.50),
  apply(terminal_prices, 2, quantile, probs = 0.99)
)
payoff_grid <- multistock_payoff_matrix(
  terminal_prices = state_grid,
  underlying = instrument_underlying,
  K = K,
  type = option_type,
  underlying_names = underlyings
)

optimization_result <- run_option_optimization_grid(
  m = NA_real_,
  v = NA_real_,
  K = K,
  type = option_type,
  option_price = option_price,
  bid_price = bid_price,
  ask_price = ask_price,
  mid_price = mid_price,
  vol_quotes = vol_quotes,
  greeks = greeks,
  S0 = S0[[1]],
  r = r,
  q = q[[1]],
  T = T,
  instrument_underlying = instrument_underlying,
  portfolio_delta_bounds = c(NA_real_, NA_real_),
  expected_payoff = horizon_inputs$expected_payoff,
  cov_payoff = horizon_inputs$cov_payoff,
  dependency_matrix = NULL,
  dependency_diagnostics = NULL,
  payoff_stats = horizon_inputs$payoff_stats,
  payoff_scenarios = horizon_inputs$payoff_scenarios,
  payoff_moment_scenarios = horizon_inputs$payoff_scenarios,
  payoff_grid = payoff_grid,
  tail_slope = rep(0.0, length(K)),
  lower_tail_slope = NULL,
  cvar_constraints = data.frame(tail_prob = numeric(), max_loss = numeric()),
  cvar_scenarios = numeric(0),
  cvar_payoff_scenarios = horizon_inputs$cvar_payoff_scenarios,
  optimization_objectives = optimization_objectives,
  risk_aversion_utility = risk_aversion_utility,
  constraint_modes = constraint_modes,
  budget = budget,
  rf_growth = exp(r * T),
  risk_aversion = risk_aversion,
  min_terminal_wealth = -Inf,
  initial_contracts = initial_contracts,
  max_invested_weight = max_invested_weight,
  force_full_investment = force_full_investment,
  constrained_optimizer = "nelder_mead",
  report_mtm_var = report_mtm_var,
  var_conf_levels = var_conf_levels,
  mtm_var_inputs_override = mtm_var_inputs,
  print_zero_weight_options = print_zero_weight_options,
  zero_weight_tol = zero_weight_tol,
  print_runs = TRUE
)

cat("Multi-stock objective summary\n")
print(optimization_result$summary_print, row.names = FALSE)

if (write_summary_csv) {
  dir <- dirname(summary_csv_file)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  write.csv(optimization_result$summary_print, summary_csv_file, row.names = FALSE)
  cat("summary_csv_file:", summary_csv_file, "\n")
}
