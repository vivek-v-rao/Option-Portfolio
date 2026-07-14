# Choose strikes for named option structures and rank the candidates.
#
# This driver is separate from xoptimize_options.r because structure selection is
# a template search over fixed option bundles, not a free portfolio optimizer.

source("option_stats.r")
options(width = 10000)

S0 <- 100
mu <- 0.03
realized_sigma <- 0.20
T <- 1.0
r <- 0.03
q <- 0.00
terminal_model <- "lognormal"
terminal_floor <- NA_real_
terminal_n_scenarios <- 401
cvar_n_scenarios <- 101

base_strikes <- c(80, 90, 100, 110, 120)
implied_sigma <- c(0.18, 0.18, 0.18, 0.18, 0.18)
structure_types <- c("straddle", "butterfly", "iron_condor")
structure_units <- 1.0
structure_objective <- "sharpe"
structure_top_n <- 10L
require_even_butterfly <- TRUE

budget <- 1000.0
risk_aversion <- 4.0
utility_gamma <- NA_real_
min_terminal_wealth <- -Inf
cvar_constraints <- data.frame(tail_prob = numeric(), max_loss = numeric())
portfolio_delta_bounds <- c(NA_real_, NA_real_)

use_bid_ask_vols <- TRUE
vol_spread_atm <- 0.01
vol_spread_quad <- 0.10
min_bid_vol <- 0.001
stock_bid_ask_spread <- if (use_bid_ask_vols) 0.01 else 0.0

report_mtm_var <- TRUE
var_horizon_days <- 1.0
var_trading_days <- 252.0
var_n_scenarios <- 501
var_return_model <- "normal_log"
var_return_mu <- mu
var_return_sigma <- realized_sigma
var_terminal_floor <- 0.0
var_conf_levels <- c(0.95, 0.99)
mtm_var_constraints <- data.frame(conf_level = numeric(), max_loss = numeric())
mtm_es_constraints <- data.frame(conf_level = numeric(), max_loss = numeric())

delta_hedge_steps <- 0L
delta_hedge_paths <- 1000L
delta_hedge_mu <- mu
delta_hedge_sigma <- realized_sigma
delta_hedge_seed <- NA_integer_
delta_hedge_stock_transaction_cost <- 0.0

write_structure_summary_csv <- FALSE
structure_summary_csv_file <- "option_structure_summary.csv"

cat("Option structure search\n")
cat("S0:", S0, "\n")
cat("mu:", mu, "\n")
cat("realized_sigma:", realized_sigma, "\n")
cat("T:", T, "\n")
cat("r:", r, "\n")
cat("q:", q, "\n")
cat("structure_types:", structure_types, "\n")
cat("structure_objective:", structure_objective, "\n")
cat("structure_units:", structure_units, "\n")
cat("structure_top_n:", structure_top_n, "\n")
cat("delta_hedge_steps:", delta_hedge_steps, "\n\n")

stopifnot(length(implied_sigma) == length(base_strikes))
K <- c(0, base_strikes, base_strikes)
option_type <- c("call", rep("call", length(base_strikes)), rep("put", length(base_strikes)))
instrument_implied_sigma <- c(NA, implied_sigma, implied_sigma)
initial_contracts <- rep(0.0, length(K))

terminal_inputs <- terminal_distribution_inputs(
  terminal_model = terminal_model,
  S0 = S0,
  mu = mu,
  sigma = realized_sigma,
  T = T,
  K = K,
  type = option_type,
  terminal_floor = terminal_floor,
  n_scenarios = terminal_n_scenarios,
  cvar_n_scenarios = cvar_n_scenarios
)

vol_quotes <- vol_bid_ask_from_quadratic_spread(
  S0 = S0,
  K = K,
  T = T,
  r = r,
  q = q,
  mid_vol = instrument_implied_sigma,
  spread_atm = if (use_bid_ask_vols) vol_spread_atm else 0.0,
  spread_quad = if (use_bid_ask_vols) vol_spread_quad else 0.0,
  min_bid_vol = min_bid_vol
)
mid_price <- bs_option_price_vec(S0, K, r, q, vol_quotes$mid_vol, T, option_type)
bid_price <- bs_option_price_vec(S0, K, r, q, vol_quotes$bid_vol, T, option_type)
ask_price <- bs_option_price_vec(S0, K, r, q, vol_quotes$ask_vol, T, option_type)
stock_idx <- which(K == 0 & option_type == "call")
if (length(stock_idx) > 0) {
  bid_price[stock_idx] <- pmax(mid_price[stock_idx] - 0.5 * stock_bid_ask_spread, 0.0)
  ask_price[stock_idx] <- mid_price[stock_idx] + 0.5 * stock_bid_ask_spread
}
option_price <- ask_price
greeks <- bs_option_greeks_table(S0, K, r, q, vol_quotes$mid_vol, T, option_type)
rf_growth <- exp(r * T)

mtm_var_inputs <- if (report_mtm_var) {
  list(
    current_mid_price = mid_price,
    S0 = S0,
    r = r,
    q = q,
    mid_vol = vol_quotes$mid_vol,
    T = T,
    var_mu = var_return_mu,
    var_sigma = var_return_sigma,
    horizon_days = var_horizon_days,
    trading_days = var_trading_days,
    n_scenarios = var_n_scenarios,
    return_model = var_return_model,
    terminal_floor = var_terminal_floor,
    conf_levels = var_conf_levels
  )
} else {
  NULL
}
delta_hedge_inputs <- if (delta_hedge_steps > 0L) {
  list(
    S0 = S0,
    r = r,
    q = q,
    mid_vol = vol_quotes$mid_vol,
    T = T,
    steps = delta_hedge_steps,
    paths = delta_hedge_paths,
    mu = delta_hedge_mu,
    sigma = delta_hedge_sigma,
    seed = delta_hedge_seed,
    stock_transaction_cost = delta_hedge_stock_transaction_cost
  )
} else {
  NULL
}

candidates <- generate_option_structure_candidates(
  K = K,
  type = option_type,
  structure_types = structure_types,
  units = structure_units,
  require_even_butterfly = require_even_butterfly
)
cat("candidate_count:", length(candidates), "\n\n")

structure_result <- evaluate_option_structures(
  candidates = candidates,
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
  T = T,
  expected_payoff = terminal_inputs$expected_payoff,
  cov_payoff = terminal_inputs$cov_payoff,
  dependency_matrix = NULL,
  payoff_stats = terminal_inputs$payoff_stats,
  payoff_scenarios = terminal_inputs$payoff_scenarios,
  payoff_moment_scenarios = NULL,
  payoff_grid = terminal_inputs$payoff_grid,
  tail_slope = terminal_inputs$tail_slope,
  cvar_constraints = cvar_constraints,
  cvar_payoff_scenarios = terminal_inputs$cvar_payoff_scenarios,
  optimization_objective = structure_objective,
  utility_gamma = utility_gamma,
  budget = budget,
  rf_growth = rf_growth,
  risk_aversion = risk_aversion,
  min_terminal_wealth = min_terminal_wealth,
  initial_contracts = initial_contracts,
  portfolio_delta_bounds = portfolio_delta_bounds,
  mtm_var_inputs = mtm_var_inputs,
  mtm_var_constraints = mtm_var_constraints,
  mtm_es_constraints = mtm_es_constraints,
  delta_hedge_inputs = delta_hedge_inputs,
  top_n = structure_top_n,
  print_zero_weight_options = FALSE
)

cat("Top option structures\n")
print(structure_result$top_summary, row.names = FALSE)

cat("\nBest structure legs\n")
print(structure_result$best_run$table_print, row.names = FALSE)

if (write_structure_summary_csv) {
  dir <- dirname(structure_summary_csv_file)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  write.csv(structure_result$summary_print, structure_summary_csv_file, row.names = FALSE)
  cat("structure_summary_csv_file:", structure_summary_csv_file, "\n")
}
