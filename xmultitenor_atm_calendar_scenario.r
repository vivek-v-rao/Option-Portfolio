# ATM two-tenor calendar-vol scenario.
#
# S0 = 100, zero rates/dividends, 1-month ATM IV = 40%, 2-month ATM IV = 30%,
# realized vol = 30%, and only the 100 strike trades at both tenors.

source("option_stats.r")
options(width = 10000)

S0 <- 100
mu <- 0.0
realized_sigma <- 0.30
r <- 0.0
q <- 0.0

expiries <- c(1.0 / 12.0, 2.0 / 12.0)
horizon <- min(expiries)
base_strikes <- c(100)
implied_sigma_by_expiry <- list(
  "0.0833333333333333" = c(0.40),
  "0.166666666666667" = c(0.30)
)

terminal_model <- "lognormal"
terminal_floor <- NA_real_
horizon_n_scenarios <- 801
cvar_n_scenarios <- 801

optimization_objectives <- c("mean_variance", "sharpe")
risk_aversion_utility <- c(1.0)
constraint_modes <- c("long_only", "nonnegative_terminal")
risk_aversion <- 4.0
budget <- 1000.0
max_invested_weight <- 1.0
force_full_investment <- TRUE
min_terminal_wealth <- 0.0
initial_contracts <- NULL

use_bid_ask_vols <- FALSE
vol_spread_atm <- 0.0
vol_spread_quad <- 0.0
min_bid_vol <- 0.001
stock_bid_ask_spread <- 0.0

print_zero_weight_options <- TRUE
zero_weight_tol <- 1e-10

cat("ATM calendar-vol multitenor scenario\n")
cat("S0:", S0, "\n")
cat("mu:", mu, "\n")
cat("realized_sigma:", realized_sigma, "\n")
cat("r:", r, "\n")
cat("q:", q, "\n")
cat("expiries:", expiries, "\n")
cat("horizon:", horizon, "\n")
cat("front_iv:", implied_sigma_by_expiry[[as.character(expiries[1])]], "\n")
cat("back_iv:", implied_sigma_by_expiry[[as.character(expiries[2])]], "\n\n")

stock_row <- data.frame(
  expiry = horizon,
  type = "call",
  strike = 0.0,
  implied_vol = NA_real_
)
option_rows <- do.call(rbind, lapply(expiries, function(expiry) {
  key <- as.character(expiry)
  vols <- implied_sigma_by_expiry[[key]]
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
mid_price <- multitenor_price_vec(S0, K, option_type, instrument_expiry, r, q, vol_quotes$mid_vol)
bid_price <- multitenor_price_vec(S0, K, option_type, instrument_expiry, r, q, vol_quotes$bid_vol)
ask_price <- multitenor_price_vec(S0, K, option_type, instrument_expiry, r, q, vol_quotes$ask_vol)
stock_idx <- which(K == 0 & option_type == "call")
if (length(stock_idx) > 0) {
  mid_price[stock_idx] <- S0
  bid_price[stock_idx] <- pmax(S0 - 0.5 * stock_bid_ask_spread, 0.0)
  ask_price[stock_idx] <- S0 + 0.5 * stock_bid_ask_spread
}
option_price <- ask_price

greeks <- multitenor_greeks_table(S0, K, option_type, instrument_expiry, r, q, vol_quotes$mid_vol)
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
    mid_price = round(mid_price, 4)
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
  min_terminal_wealth = min_terminal_wealth,
  initial_contracts = initial_contracts,
  max_invested_weight = max_invested_weight,
  force_full_investment = force_full_investment,
  constrained_optimizer = "nelder_mead",
  constrained_optimizer_max_starts = 10,
  constrained_optimizer_max_iter = 1000,
  report_mtm_var = FALSE,
  print_zero_weight_options = print_zero_weight_options,
  zero_weight_tol = zero_weight_tol,
  print_runs = TRUE
)

cat("ATM calendar-vol objective summary\n")
print(optimization_result$summary_print, row.names = FALSE)

cat("\nInterpretation hint:\n")
cat("Under nonnegative_terminal constraints, a positive 2-month straddle position\n")
cat("paired with a negative 1-month straddle position is the expected calendar-vol trade.\n")
