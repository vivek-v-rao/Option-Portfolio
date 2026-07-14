# Multi-tenor option portfolio prototype.
#
# The optimizer horizon is the first expiry. Options expiring at that horizon
# settle to payoff; later-tenor options are marked at the horizon with
# Black-Scholes using remaining maturity and sticky-strike implied vols.

source("option_stats.r")
options(width = 10000)

S0 <- 100
mu <- 0.03
realized_sigma <- 0.20
r <- 0.03
q <- 0.00

expiries <- c(0.25, 0.50)
horizon <- min(expiries)
base_strikes <- c(90, 100, 110)
implied_sigma_by_expiry <- list(
  "0.25" = c(0.19, 0.18, 0.19),
  "0.5" = c(0.20, 0.19, 0.20)
)

terminal_model <- "lognormal"
terminal_floor <- NA_real_
horizon_n_scenarios <- 121
cvar_n_scenarios <- 121

optimization_objectives <- c("mean_variance")
risk_aversion_utility <- c(1.0)
constraint_modes <- c("long_only")
risk_aversion <- 4.0
budget <- 1000.0
max_invested_weight <- 1.0
force_full_investment <- TRUE
initial_contracts <- NULL

use_bid_ask_vols <- TRUE
vol_spread_atm <- 0.01
vol_spread_quad <- 0.10
min_bid_vol <- 0.001
stock_bid_ask_spread <- if (use_bid_ask_vols) 0.01 else 0.0

print_zero_weight_options <- FALSE
zero_weight_tol <- 1e-10
write_summary_csv <- FALSE
summary_csv_file <- "multitenor_objective_summary.csv"

cat("Multi-tenor option optimization\n")
cat("S0:", S0, "\n")
cat("mu:", mu, "\n")
cat("realized_sigma:", realized_sigma, "\n")
cat("expiries:", expiries, "\n")
cat("horizon:", horizon, "\n")
cat("r:", r, "\n")
cat("q:", q, "\n")
cat("optimization_objectives:", optimization_objectives, "\n")
cat("constraint_modes:", constraint_modes, "\n\n")

stock_row <- data.frame(
  expiry = horizon,
  type = "call",
  strike = 0.0,
  implied_vol = NA_real_
)
option_rows <- do.call(rbind, lapply(expiries, function(expiry) {
  key <- as.character(expiry)
  vols <- implied_sigma_by_expiry[[key]]
  if (is.null(vols)) {
    stop("Missing implied vols for expiry ", key)
  }
  if (length(vols) != length(base_strikes)) {
    stop("implied vols for expiry ", key, " must match base_strikes length")
  }
  rbind(
    data.frame(expiry = expiry, type = "call", strike = base_strikes, implied_vol = vols),
    data.frame(expiry = expiry, type = "put", strike = base_strikes, implied_vol = vols)
  )
}))
instrument_table <- rbind(stock_row, option_rows)
K <- instrument_table$strike
option_type <- instrument_table$type
instrument_expiry <- instrument_table$expiry
instrument_implied_sigma <- instrument_table$implied_vol

if (is.null(initial_contracts)) {
  initial_contracts <- rep(0.0, nrow(instrument_table))
}

terminal_inputs <- terminal_distribution_inputs(
  terminal_model = terminal_model,
  S0 = S0,
  mu = mu,
  sigma = realized_sigma,
  T = horizon,
  K = K,
  type = option_type,
  terminal_floor = terminal_floor,
  n_scenarios = horizon_n_scenarios,
  cvar_n_scenarios = cvar_n_scenarios
)

vol_quotes <- data.frame(
  mid_vol = instrument_implied_sigma,
  bid_vol = instrument_implied_sigma,
  ask_vol = instrument_implied_sigma,
  vol_spread = rep(0.0, length(K))
)
for (expiry in expiries) {
  idx <- which(instrument_expiry == expiry & K > 0)
  one <- vol_bid_ask_from_quadratic_spread(
    S0 = S0,
    K = K[idx],
    T = expiry,
    r = r,
    q = q,
    mid_vol = instrument_implied_sigma[idx],
    spread_atm = if (use_bid_ask_vols) vol_spread_atm else 0.0,
    spread_quad = if (use_bid_ask_vols) vol_spread_quad else 0.0,
    min_bid_vol = min_bid_vol
  )
  vol_quotes[idx, ] <- one
}

mid_price <- multitenor_price_vec(
  S0 = S0,
  K = K,
  type = option_type,
  expiry = instrument_expiry,
  r = r,
  q = q,
  sigma = vol_quotes$mid_vol
)
bid_price <- multitenor_price_vec(S0, K, option_type, instrument_expiry, r, q, vol_quotes$bid_vol)
ask_price <- multitenor_price_vec(S0, K, option_type, instrument_expiry, r, q, vol_quotes$ask_vol)
stock_idx <- which(K == 0 & option_type == "call")
if (length(stock_idx) > 0) {
  mid_price[stock_idx] <- S0
  bid_price[stock_idx] <- pmax(S0 - 0.5 * stock_bid_ask_spread, 0.0)
  ask_price[stock_idx] <- S0 + 0.5 * stock_bid_ask_spread
}
option_price <- ask_price

greeks <- multitenor_greeks_table(
  S0 = S0,
  K = K,
  type = option_type,
  expiry = instrument_expiry,
  r = r,
  q = q,
  sigma = vol_quotes$mid_vol
)
if (length(stock_idx) > 0) {
  greeks$delta[stock_idx] <- 1.0
  greeks$gamma[stock_idx] <- 0.0
  greeks$vega[stock_idx] <- 0.0
}

horizon_inputs <- multitenor_horizon_inputs(
  horizon_prices = terminal_inputs$scenario_prices,
  cvar_prices = terminal_inputs$cvar_scenarios,
  state_grid = terminal_inputs$state_grid,
  K = K,
  type = option_type,
  expiry = instrument_expiry,
  horizon = horizon,
  r = r,
  q = q,
  sigma = vol_quotes$mid_vol
)

cat("Instrument universe\n")
print(
  data.frame(
    id = seq_along(K),
    expiry = instrument_expiry,
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

optimization_result <- run_option_optimization_grid(
  m = terminal_inputs$m,
  v = terminal_inputs$v,
  K = K,
  type = option_type,
  option_price = option_price,
  bid_price = bid_price,
  ask_price = ask_price,
  mid_price = mid_price,
  vol_quotes = vol_quotes,
  greeks = greeks,
  S0 = S0,
  r = r,
  q = q,
  T = horizon,
  instrument_expiry = instrument_expiry,
  portfolio_delta_bounds = c(NA_real_, NA_real_),
  expected_payoff = horizon_inputs$expected_value,
  cov_payoff = horizon_inputs$cov_value,
  dependency_matrix = NULL,
  dependency_diagnostics = NULL,
  payoff_stats = horizon_inputs$value_stats,
  payoff_scenarios = horizon_inputs$value_scenarios,
  payoff_moment_scenarios = horizon_inputs$value_scenarios,
  payoff_grid = horizon_inputs$value_grid,
  tail_slope = option_tail_slope(K, option_type),
  lower_tail_slope = NULL,
  cvar_constraints = data.frame(tail_prob = numeric(), max_loss = numeric()),
  cvar_scenarios = terminal_inputs$cvar_scenarios,
  cvar_payoff_scenarios = horizon_inputs$cvar_value_scenarios,
  optimization_objectives = optimization_objectives,
  risk_aversion_utility = risk_aversion_utility,
  constraint_modes = constraint_modes,
  budget = budget,
  rf_growth = exp(r * horizon),
  risk_aversion = risk_aversion,
  min_terminal_wealth = -Inf,
  initial_contracts = initial_contracts,
  max_invested_weight = max_invested_weight,
  force_full_investment = force_full_investment,
  constrained_optimizer = "nelder_mead",
  report_mtm_var = FALSE,
  print_zero_weight_options = print_zero_weight_options,
  zero_weight_tol = zero_weight_tol,
  print_runs = TRUE
)

cat("Multi-tenor objective summary\n")
print(optimization_result$summary_print, row.names = FALSE)

if (write_summary_csv) {
  dir <- dirname(summary_csv_file)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  write.csv(optimization_result$summary_print, summary_csv_file, row.names = FALSE)
  cat("summary_csv_file:", summary_csv_file, "\n")
}
