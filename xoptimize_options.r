# Optimize a single-underlying option portfolio under multiple objectives and
# constraint modes.
#
# The instrument universe contains the stock, reported as a call with strike 0,
# plus calls and puts at base_strikes. Calls and puts use the same strike-level
# implied volatility vector. Payoff moments are computed under a configurable
# terminal stock distribution with user-specified realized volatility. Price-
# space distributions can use terminal_floor to collapse mass below the floor
# into a point mass.
#
# For each constraint mode and objective, the script prints an option table with
# executable weights/contracts, single-instrument payoff statistics, and
# portfolio correlation. It also prints portfolio wealth statistics: mean, sd,
# Sharpe, skew, excess kurtosis, minimum terminal wealth, tail slope, and ES
# losses at the configured CVaR tail probabilities.
#
# Supported objectives include Sharpe, mean-variance, adjusted Sharpe,
# dependency-penalty, and expected utility over risk_aversion_utility values.
# The dependency-penalty objective uses a scenario-estimated conditional-payout
# dependency matrix with eigenvalue flooring before applying the quadratic
# penalty. Supported constraint modes include long-only and nonnegative terminal
# wealth; the latter also enforces min_terminal_wealth and configured CVaR/ES
# loss caps. When the SLSQP path is available, ES caps are native nonlinear
# inequality constraints; the base-R fallback retains the existing ES penalty
# treatment.
#
# Option transaction costs can be modeled with bid/ask implied vols. The default
# spread model widens vol spreads quadratically in log-moneyness; the stock row
# has a separate bid/ask spread. If initial_contracts is nonzero, transaction
# costs are charged only on trades away from the initial portfolio.
#
# For nonnegative-terminal nonlinear objectives, constrained_optimizer controls
# the constrained optimizer. "auto" uses nloptr SLSQP with analytic C++
# gradients when available and otherwise falls back to the base-R Nelder-Mead
# constrained search. "hybrid" runs both and keeps the better feasible result.
# The summary table reports optimizer, status code, evaluation count, and whether
# fallback was used.
# constrained_optimizer_use_warm_starts can pass prior nonnegative-terminal
# solutions into later objectives; it is off by default because it is scenario-
# dependent and can increase SLSQP work.
#
# The final Optimization objective summary and per-run raw option tables can
# optionally be written to CSV via write_summary_csv, write_portfolio_tables_csv,
# and their file/directory settings.

source("option_stats.r")
options(width = 10000)
start_time <- proc.time()
section_timings <- data.frame(section = character(), elapsed_seconds = numeric())

time_section <- function(section, expr) {
  section_start <- proc.time()
  value <- force(expr)
  section_elapsed <- proc.time() - section_start
  section_timings <<- rbind(
    section_timings,
    data.frame(section = section, elapsed_seconds = as.numeric(section_elapsed["elapsed"]))
  )
  value
}

append_section_timings <- function(prefix, timings) {
  if (!is.null(timings) && nrow(timings) > 0) {
    section_timings <<- rbind(
      section_timings,
      data.frame(
        section = paste0(prefix, ":", timings$section),
        elapsed_seconds = timings$elapsed_seconds
      )
    )
  }
}

ensure_parent_dir <- function(file) {
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  invisible(file)
}

parse_option_cli_args <- function(args) {
  out <- list(config_file = NA_character_, scenario = NA_character_)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--config") {
      i <- i + 1L
      if (i > length(args)) stop("--config requires a file path")
      out$config_file <- args[[i]]
    } else if (startsWith(arg, "--config=")) {
      out$config_file <- sub("^--config=", "", arg)
    } else if (arg == "--scenario") {
      i <- i + 1L
      if (i > length(args)) stop("--scenario requires a preset name")
      out$scenario <- args[[i]]
    } else if (startsWith(arg, "--scenario=")) {
      out$scenario <- sub("^--scenario=", "", arg)
    } else {
      stop("Unknown command-line argument: ", arg)
    }
    i <- i + 1L
  }
  out
}

coerce_config_overrides <- function(overrides) {
  null_as_na <- c(
    "terminal_floor",
    "min_terminal_wealth",
    "stock_bid_ask_spread",
    "plot_constraint_mode",
    "plot_objective",
    "plot_utility_gamma",
    "portfolio_delta_bounds",
    "var_return_mu",
    "var_return_sigma",
    "var_terminal_floor",
    "max_active_option_positions",
    "delta_hedge_mu",
    "delta_hedge_sigma",
    "delta_hedge_seed"
  )
  for (nm in names(overrides)) {
    if (is.null(overrides[[nm]]) && nm %in% null_as_na) {
      overrides[[nm]] <- if (nm == "portfolio_delta_bounds") c(NA_real_, NA_real_) else NA_real_
    }
  }
  if ("cvar_constraints" %in% names(overrides)) {
    x <- overrides$cvar_constraints
    if (is.null(x)) {
      overrides$cvar_constraints <- data.frame(tail_prob = numeric(), max_loss = numeric())
    } else if (is.data.frame(x)) {
      overrides$cvar_constraints <- data.frame(
        tail_prob = as.numeric(x$tail_prob),
        max_loss = as.numeric(x$max_loss)
      )
    } else if (is.list(x) && all(c("tail_prob", "max_loss") %in% names(x))) {
      overrides$cvar_constraints <- data.frame(
        tail_prob = as.numeric(x$tail_prob),
        max_loss = as.numeric(x$max_loss)
      )
    } else {
      stop("cvar_constraints in config file must contain tail_prob and max_loss")
    }
  }
  for (constraint_name in c("mtm_var_constraints", "mtm_es_constraints")) {
    if (constraint_name %in% names(overrides)) {
      x <- overrides[[constraint_name]]
      if (is.null(x)) {
        overrides[[constraint_name]] <- data.frame(conf_level = numeric(), max_loss = numeric())
      } else if (is.data.frame(x)) {
        overrides[[constraint_name]] <- data.frame(
          conf_level = as.numeric(x$conf_level),
          max_loss = as.numeric(x$max_loss)
        )
      } else if (is.list(x) && all(c("conf_level", "max_loss") %in% names(x))) {
        overrides[[constraint_name]] <- data.frame(
          conf_level = as.numeric(x$conf_level),
          max_loss = as.numeric(x$max_loss)
        )
      } else {
        stop(constraint_name, " in config file must contain conf_level and max_loss")
      }
    }
  }
  overrides
}

read_option_config_file <- function(file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Reading JSON config files requires the jsonlite package")
  }
  raw <- jsonlite::fromJSON(file, simplifyVector = TRUE)
  scenario <- if ("scenario" %in% names(raw)) raw$scenario else NA_character_
  overrides <- if ("config" %in% names(raw)) raw$config else raw[setdiff(names(raw), "scenario")]
  if (is.null(overrides)) {
    overrides <- list()
  }
  list(
    scenario = scenario,
    overrides = coerce_config_overrides(overrides)
  )
}

cli_args <- parse_option_cli_args(commandArgs(trailingOnly = TRUE))
option_config_file_path <- if (!is.na(cli_args$config_file)) {
  cli_args$config_file
} else if (exists("option_config_file", inherits = FALSE)) {
  option_config_file
} else {
  NA_character_
}
file_config <- if (!is.na(option_config_file_path)) {
  read_option_config_file(option_config_file_path)
} else {
  list(scenario = NA_character_, overrides = list())
}

scenario <- "full"
if (!is.na(file_config$scenario)) scenario <- file_config$scenario
if (exists("option_scenario", inherits = FALSE)) scenario <- option_scenario
if (!is.na(cli_args$scenario)) scenario <- cli_args$scenario

config <- list(
  S0 = 100,
  mu = 0.03,
  realized_sigma = 0.20,
  T = 1.0,
  r = 0.03,
  q = 0.00,
  terminal_model = "lognormal",
  terminal_floor = NA_real_,
  terminal_n_scenarios = 401,
  mixture_weights = c(0.85, 0.15),
  mixture_mu = c(0.06, -0.25),
  mixture_sigma = c(0.16, 0.45),
  hyperbolic_alpha = 10.0,
  hyperbolic_beta = -2.0,
  gh_lambda = 1.0,
  gh_alpha = 10.0,
  gh_beta = -2.0,
  nig_alpha = 10.0,
  nig_beta = -2.0,
  vg_shape = 1.0,
  vg_beta = -0.2,
  risk_aversion = 4.0,
  budget = 1000.0,
  optimization_objectives = c("sharpe", "mean_variance", "dependency_penalty", "expected_utility"),
  risk_aversion_utility = c(0.0, 0.001, 1.0, 2.0),
  min_terminal_wealth = NA_real_,
  cvar_constraints = data.frame(
    tail_prob = c(0.01, 0.05),
    max_loss = c(200.0, 100.0)
  ),
  cvar_n_scenarios = 401,
  constraint_modes = c("long_only", "nonnegative_terminal"),
  use_bid_ask_vols = TRUE,
  vol_spread_atm = 0.01,
  vol_spread_quad = 0.10,
  min_bid_vol = 0.001,
  stock_bid_ask_spread = NA_real_,
  max_invested_weight = 1.0,
  portfolio_delta_bounds = c(NA_real_, NA_real_),
  force_full_investment = TRUE, # FALSE
  max_iter = 2000,
  nonnegative_expected_utility_max_starts = Inf,
  nonnegative_expected_utility_max_iter = 2000,
  constrained_optimizer_max_starts = Inf,
  constrained_optimizer_max_iter = 2000,
  constrained_optimizer_use_warm_starts = FALSE,
  constrained_optimizer = "auto",
  report_mtm_var = TRUE,
  var_horizon_days = 1.0,
  var_trading_days = 252.0,
  var_n_scenarios = 1001,
  var_return_model = "normal_log",
  var_return_mu = NA_real_,
  var_return_sigma = NA_real_,
  var_terminal_floor = 0.0,
  var_conf_levels = c(0.95, 0.99),
  mtm_var_constraints = data.frame(conf_level = numeric(), max_loss = numeric()),
  mtm_es_constraints = data.frame(conf_level = numeric(), max_loss = numeric()),
  delta_hedge_steps = 0L,
  delta_hedge_paths = 1000L,
  delta_hedge_mu = NA_real_,
  delta_hedge_sigma = NA_real_,
  delta_hedge_seed = NA_integer_,
  delta_hedge_stock_transaction_cost = 0.0,
  max_active_option_positions = NA_integer_,
  prune_positions_by = "position_weight",
  prune_repair_constraints = TRUE,
  integer_contracts = FALSE,
  integer_rounding_neighborhood = 1L,
  integer_max_search_instruments = 6L,
  print_zero_weight_options = TRUE,
  zero_weight_tol = 1e-12,
  report_greeks = TRUE,
  run_rebalance_sanity_checks = TRUE,
  write_summary_csv = FALSE,
  summary_csv_file = "optimization_objective_summary.csv",
  write_portfolio_tables_csv = FALSE,
  portfolio_tables_csv_dir = "portfolio_tables",
  write_combined_portfolio_table_csv = FALSE,
  combined_portfolio_table_csv_file = "portfolio_tables/portfolio_tables_combined.csv",
  write_portfolio_plot = FALSE,
  portfolio_plot_file = "option_portfolio_report.png",
  plot_constraint_mode = NA_character_,
  plot_objective = NA_character_,
  plot_utility_gamma = NA_real_,
  base_strikes = c(80, 90, 100, 110, 120),
  implied_sigma = c(0.18, 0.18, 0.18, 0.18, 0.18),
  initial_contracts = NULL,
  tol = 1e-10
)

apply_config_overrides <- function(config, overrides, label) {
  if ("constrained_expected_utility_optimizer" %in% names(overrides) &&
      !"constrained_optimizer" %in% names(overrides)) {
    overrides$constrained_optimizer <- overrides$constrained_expected_utility_optimizer
    overrides$constrained_expected_utility_optimizer <- NULL
  }
  unknown_names <- setdiff(names(overrides), names(config))
  if (length(unknown_names) > 0) {
    stop("Unknown ", label, " names: ", paste(unknown_names, collapse = ", "))
  }
  for (nm in names(overrides)) {
    config[[nm]] <- overrides[[nm]]
  }
  config
}

config_presets <- list(
  full = list(),
  quick = list(
    optimization_objectives = c("sharpe", "mean_variance", "dependency_penalty", "expected_utility"),
    risk_aversion_utility = c(1.0),
    cvar_n_scenarios = 81,
    print_zero_weight_options = FALSE,
    zero_weight_tol = 1e-8,
    write_summary_csv = FALSE,
    write_portfolio_tables_csv = FALSE,
    write_combined_portfolio_table_csv = FALSE,
    write_portfolio_plot = FALSE,
    base_strikes = c(90, 100, 110),
    implied_sigma = c(0.18, 0.18, 0.18)
  ),
  no_costs = list(
    use_bid_ask_vols = FALSE,
    stock_bid_ask_spread = NA_real_
  ),
  long_only = list(
    constraint_modes = c("long_only")
  ),
  long_short = list(
    constraint_modes = c("long_short"),
    optimization_objectives = c("mean_variance", "sharpe")
  ),
  logistic = list(
    terminal_model = "logistic",
    terminal_floor = NA_real_
  ),
  hyperbolic_secant = list(
    terminal_model = "hyperbolic_secant",
    terminal_floor = NA_real_
  ),
  symmetric_hyperbolic = list(
    terminal_model = "symmetric_hyperbolic",
    terminal_floor = NA_real_
  ),
  hyperbolic = list(
    terminal_model = "hyperbolic",
    terminal_floor = NA_real_,
    hyperbolic_beta = -2.0
  ),
  generalized_hyperbolic = list(
    terminal_model = "generalized_hyperbolic",
    terminal_floor = NA_real_
  ),
  normal_inverse_gaussian = list(
    terminal_model = "normal_inverse_gaussian",
    terminal_floor = NA_real_
  ),
  variance_gamma = list(
    terminal_model = "variance_gamma",
    terminal_floor = NA_real_
  ),
  normal_floor = list(
    terminal_model = "normal",
    terminal_floor = 0.0
  ),
  normal_unbounded = list(
    terminal_model = "normal",
    terminal_floor = NA_real_
  ),
  crash_mixture = list(
    terminal_model = "lognormal_mixture",
    terminal_floor = NA_real_,
    terminal_n_scenarios = 801,
    cvar_n_scenarios = 801,
    nonnegative_expected_utility_max_starts = 5,
    nonnegative_expected_utility_max_iter = 1000,
    mixture_weights = c(0.85, 0.15),
    mixture_mu = c(0.06, -0.25),
    mixture_sigma = c(0.16, 0.45),
    implied_sigma = rep(0.18, length(config$base_strikes))
  ),
  tail_risk = list(
    constraint_modes = c("nonnegative_terminal"),
    optimization_objectives = c("mean_variance", "expected_utility"),
    risk_aversion_utility = c(1.0, 2.0),
    cvar_constraints = data.frame(
      tail_prob = c(0.01, 0.05),
      max_loss = c(100.0, 50.0)
    )
  )
)

if (length(scenario) != 1 || is.na(scenario) || !scenario %in% names(config_presets)) {
  stop("Unknown option scenario: ", paste(scenario, collapse = ", "),
       ". Valid scenarios: ", paste(names(config_presets), collapse = ", "))
}
config <- apply_config_overrides(config, config_presets[[scenario]], "config preset")

if (length(file_config$overrides) > 0) {
  config <- apply_config_overrides(config, file_config$overrides, "config file")
}

if (exists("option_config_overrides", inherits = FALSE)) {
  config <- apply_config_overrides(
    config,
    coerce_config_overrides(option_config_overrides),
    "option_config_overrides"
  )
}

if (is.na(config$min_terminal_wealth)) {
  config$min_terminal_wealth <- if (scenario == "tail_risk") 0.90 * config$budget else 1e-6 * config$budget
}
if (is.na(config$stock_bid_ask_spread)) {
  config$stock_bid_ask_spread <- if (config$use_bid_ask_vols) 0.01 else 0.0
}
if (is.na(config$var_return_mu)) {
  config$var_return_mu <- config$mu
}
if (is.na(config$var_return_sigma)) {
  config$var_return_sigma <- config$realized_sigma
}
if (is.na(config$delta_hedge_mu)) {
  config$delta_hedge_mu <- config$mu
}
if (is.na(config$delta_hedge_sigma)) {
  config$delta_hedge_sigma <- config$realized_sigma
}
config$scenario <- scenario
invisible(list2env(config, envir = environment()))

cat("scenario:", scenario,
	"\nconfig_file:", option_config_file_path,
	"\nS0:", S0, "\nmu:", mu, "\nrealized_vol", realized_sigma, "\nT:", T, "\nr:", r,
	"\nq:", q, "\nrisk_aversion:", risk_aversion, "\nbudget:", budget,
	"\ncpp_kernels_loaded:", option_kernels_loaded,
	"\nterminal_model:", terminal_model,
	"\nterminal_floor:", terminal_floor,
	"\nterminal_n_scenarios:", terminal_n_scenarios,
	"\nmixture_weights:", mixture_weights,
	"\nmixture_mu:", mixture_mu,
	"\nmixture_sigma:", mixture_sigma,
	"\nhyperbolic_alpha:", hyperbolic_alpha,
	"\nhyperbolic_beta:", hyperbolic_beta,
	"\ngh_lambda:", gh_lambda,
	"\ngh_alpha:", gh_alpha,
	"\ngh_beta:", gh_beta,
	"\nnig_alpha:", nig_alpha,
	"\nnig_beta:", nig_beta,
	"\nvg_shape:", vg_shape,
	"\nvg_beta:", vg_beta,
	"\noptimization_objectives:", optimization_objectives,
	"\nrisk_aversion_utility:", risk_aversion_utility,
	"\nmin_terminal_wealth:", min_terminal_wealth,
	"\ncvar_n_scenarios:", cvar_n_scenarios,
	"\nconstraint_modes:", constraint_modes,
	"\nuse_bid_ask_vols:", use_bid_ask_vols,
	"\nvol_spread_atm:", vol_spread_atm,
	"\nvol_spread_quad:", vol_spread_quad,
	"\nmin_bid_vol:", min_bid_vol,
	"\nstock_bid_ask_spread:", stock_bid_ask_spread,
	"\nmax_invested_weight:", max_invested_weight,
	"\nportfolio_delta_bounds:", portfolio_delta_bounds,
	"\nforce_full_investment:", force_full_investment, "\nmax_iter:",
	max_iter,
	"\nnonnegative_expected_utility_max_starts:", nonnegative_expected_utility_max_starts,
	"\nnonnegative_expected_utility_max_iter:", nonnegative_expected_utility_max_iter,
	"\nconstrained_optimizer_max_starts:", constrained_optimizer_max_starts,
	"\nconstrained_optimizer_max_iter:", constrained_optimizer_max_iter,
	"\nconstrained_optimizer_use_warm_starts:", constrained_optimizer_use_warm_starts,
	"\nconstrained_optimizer:", constrained_optimizer,
	"\nreport_mtm_var:", report_mtm_var,
	"\nvar_horizon_days:", var_horizon_days,
	"\nvar_trading_days:", var_trading_days,
	"\nvar_n_scenarios:", var_n_scenarios,
	"\nvar_return_model:", var_return_model,
	"\nvar_return_mu:", var_return_mu,
	"\nvar_return_sigma:", var_return_sigma,
	"\nvar_terminal_floor:", var_terminal_floor,
	"\nvar_conf_levels:", var_conf_levels,
	"\nmtm_var_constraints_n:", nrow(mtm_var_constraints),
	"\nmtm_es_constraints_n:", nrow(mtm_es_constraints),
	"\ndelta_hedge_steps:", delta_hedge_steps,
	"\ndelta_hedge_paths:", delta_hedge_paths,
	"\ndelta_hedge_mu:", delta_hedge_mu,
	"\ndelta_hedge_sigma:", delta_hedge_sigma,
	"\ndelta_hedge_seed:", delta_hedge_seed,
	"\ndelta_hedge_stock_transaction_cost:", delta_hedge_stock_transaction_cost,
	"\nmax_active_option_positions:", max_active_option_positions,
	"\nprune_positions_by:", prune_positions_by,
	"\nprune_repair_constraints:", prune_repair_constraints,
	"\ninteger_contracts:", integer_contracts,
	"\ninteger_rounding_neighborhood:", integer_rounding_neighborhood,
	"\ninteger_max_search_instruments:", integer_max_search_instruments,
	"\nprint_zero_weight_options:", print_zero_weight_options,
	"\nzero_weight_tol:", zero_weight_tol,
	"\nreport_greeks:", report_greeks,
	"\nrun_rebalance_sanity_checks:", run_rebalance_sanity_checks,
	"\nwrite_summary_csv:", write_summary_csv,
	"\nsummary_csv_file:", summary_csv_file,
	"\nwrite_portfolio_tables_csv:", write_portfolio_tables_csv,
	"\nportfolio_tables_csv_dir:", portfolio_tables_csv_dir,
	"\nwrite_combined_portfolio_table_csv:", write_combined_portfolio_table_csv,
	"\ncombined_portfolio_table_csv_file:", combined_portfolio_table_csv_file,
	"\nwrite_portfolio_plot:", write_portfolio_plot,
	"\nportfolio_plot_file:", portfolio_plot_file,
	"\nplot_constraint_mode:", plot_constraint_mode,
	"\nplot_objective:", plot_objective,
	"\nplot_utility_gamma:", plot_utility_gamma, "\n\n")
if (run_rebalance_sanity_checks) {
  invisible(time_section("rebalance_sanity_checks", check_rebalance_trade_costs()))
}
cat("cvar_constraints:\n")
print(cvar_constraints, row.names = FALSE)
cat("\n")
cat("mtm_var_constraints:\n")
print(mtm_var_constraints, row.names = FALSE)
cat("\n")
cat("mtm_es_constraints:\n")
print(mtm_es_constraints, row.names = FALSE)
cat("\n")
stopifnot(length(implied_sigma) == length(base_strikes))
K <- c(0, base_strikes, base_strikes)
option_type <- c("call", rep("call", length(base_strikes)), rep("put", length(base_strikes)))
instrument_implied_sigma <- c(NA, implied_sigma, implied_sigma)
if (is.null(initial_contracts)) {
  initial_contracts <- rep(0.0, length(K))
} else {
  initial_contracts <- as.numeric(initial_contracts)
  if (length(initial_contracts) != length(K)) {
    stop("initial_contracts must have length ", length(K),
         " for stock + calls + puts at configured strikes")
  }
}

terminal_inputs <- time_section("terminal_distribution_inputs", terminal_distribution_inputs(
  terminal_model = terminal_model,
  S0 = S0,
  mu = mu,
  sigma = realized_sigma,
  T = T,
  K = K,
  type = option_type,
  terminal_floor = terminal_floor,
  n_scenarios = terminal_n_scenarios,
  cvar_n_scenarios = cvar_n_scenarios,
  mixture_weights = mixture_weights,
  mixture_mu = mixture_mu,
  mixture_sigma = mixture_sigma,
  hyperbolic_alpha = hyperbolic_alpha,
  hyperbolic_beta = hyperbolic_beta,
  gh_lambda = gh_lambda,
  gh_alpha = gh_alpha,
  gh_beta = gh_beta,
  nig_alpha = nig_alpha,
  nig_beta = nig_beta,
  vg_shape = vg_shape,
  vg_beta = vg_beta
))
append_section_timings("terminal_distribution_inputs", terminal_inputs$section_timings)
m <- terminal_inputs$m
v <- terminal_inputs$v
payoff_scenarios <- terminal_inputs$payoff_scenarios
payoff_moment_scenarios <- if (terminal_model == "lognormal" && is.na(terminal_floor)) NULL else payoff_scenarios
payoff_stats <- terminal_inputs$payoff_stats
expected_payoff <- terminal_inputs$expected_payoff
cov_payoff <- terminal_inputs$cov_payoff
payoff_grid <- terminal_inputs$payoff_grid
tail_slope <- terminal_inputs$tail_slope
lower_tail_slope <- terminal_inputs$lower_tail_slope
cvar_scenarios <- terminal_inputs$cvar_scenarios
cvar_payoff_scenarios <- terminal_inputs$cvar_payoff_scenarios
dependency_inputs <- time_section("payoff_dependency_matrix", {
  dependency_raw <- option_payoff_dependency_matrix(payoff_scenarios)
  psd_repair_matrix(dependency_raw, eigen_floor = 1e-8)
})
dependency_matrix <- dependency_inputs$matrix
dependency_diagnostics <- dependency_inputs
market_inputs <- time_section("market_prices_and_vols", {
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
mid_price <- bs_option_price_vec(
  S0,
  K,
  r,
  q,
  vol_quotes$mid_vol,
  T,
  option_type
)
bid_price <- bs_option_price_vec(S0, K, r, q, vol_quotes$bid_vol, T, option_type)
ask_price <- bs_option_price_vec(S0, K, r, q, vol_quotes$ask_vol, T, option_type)
stock_idx <- which(K == 0 & option_type == "call")
if (length(stock_idx) > 0) {
  bid_price[stock_idx] <- pmax(mid_price[stock_idx] - 0.5 * stock_bid_ask_spread, 0.0)
  ask_price[stock_idx] <- mid_price[stock_idx] + 0.5 * stock_bid_ask_spread
}
option_price <- ask_price
greeks <- bs_option_greeks_table(S0, K, r, q, vol_quotes$mid_vol, T, option_type)
list(
  vol_quotes = vol_quotes,
  mid_price = mid_price,
  bid_price = bid_price,
  ask_price = ask_price,
  option_price = option_price,
  greeks = greeks
)
})
vol_quotes <- market_inputs$vol_quotes
mid_price <- market_inputs$mid_price
bid_price <- market_inputs$bid_price
ask_price <- market_inputs$ask_price
option_price <- market_inputs$option_price
greeks <- market_inputs$greeks
cat("initial_contracts:\n")
print(
  data.frame(
    type = option_type,
    strike = K,
    initial_contracts = initial_contracts
  ),
  row.names = FALSE
)
cat("\n")
rf_growth <- exp(r * T)

optimization_result <- time_section("optimization_grid", run_option_optimization_grid(
  m = m,
  v = v,
  K = K,
  type = option_type,
  option_price = option_price,
  bid_price = bid_price,
  ask_price = ask_price,
  mid_price = mid_price,
  vol_quotes = vol_quotes,
  greeks = if (report_greeks) greeks else NULL,
  S0 = S0,
  r = r,
  q = q,
  T = T,
  expected_payoff = expected_payoff,
  cov_payoff = cov_payoff,
  dependency_matrix = dependency_matrix,
  dependency_diagnostics = dependency_diagnostics,
  payoff_stats = payoff_stats,
  payoff_scenarios = payoff_scenarios,
  payoff_moment_scenarios = payoff_moment_scenarios,
  payoff_grid = payoff_grid,
  tail_slope = tail_slope,
  lower_tail_slope = lower_tail_slope,
  cvar_constraints = cvar_constraints,
  cvar_scenarios = cvar_scenarios,
  cvar_payoff_scenarios = cvar_payoff_scenarios,
  optimization_objectives = optimization_objectives,
  risk_aversion_utility = risk_aversion_utility,
  constraint_modes = constraint_modes,
  budget = budget,
  rf_growth = rf_growth,
  risk_aversion = risk_aversion,
  min_terminal_wealth = min_terminal_wealth,
  initial_contracts = initial_contracts,
  max_invested_weight = max_invested_weight,
  portfolio_delta_bounds = portfolio_delta_bounds,
  force_full_investment = force_full_investment,
  nonnegative_expected_utility_max_starts = nonnegative_expected_utility_max_starts,
  nonnegative_expected_utility_max_iter = nonnegative_expected_utility_max_iter,
  constrained_optimizer_max_starts = constrained_optimizer_max_starts,
  constrained_optimizer_max_iter = constrained_optimizer_max_iter,
  constrained_optimizer_use_warm_starts = constrained_optimizer_use_warm_starts,
  constrained_optimizer = constrained_optimizer,
  report_mtm_var = report_mtm_var,
  var_horizon_days = var_horizon_days,
  var_trading_days = var_trading_days,
  var_n_scenarios = var_n_scenarios,
  var_return_model = var_return_model,
  var_return_mu = var_return_mu,
  var_return_sigma = var_return_sigma,
  var_terminal_floor = var_terminal_floor,
  var_conf_levels = var_conf_levels,
  mtm_var_constraints = mtm_var_constraints,
  mtm_es_constraints = mtm_es_constraints,
  delta_hedge_steps = delta_hedge_steps,
  delta_hedge_paths = delta_hedge_paths,
  delta_hedge_mu = delta_hedge_mu,
  delta_hedge_sigma = delta_hedge_sigma,
  delta_hedge_seed = delta_hedge_seed,
  delta_hedge_stock_transaction_cost = delta_hedge_stock_transaction_cost,
  max_active_option_positions = max_active_option_positions,
  prune_positions_by = prune_positions_by,
  prune_repair_constraints = prune_repair_constraints,
  integer_contracts = integer_contracts,
  integer_rounding_neighborhood = integer_rounding_neighborhood,
  integer_max_search_instruments = integer_max_search_instruments,
  print_zero_weight_options = print_zero_weight_options,
  zero_weight_tol = zero_weight_tol,
  tol = tol,
  print_runs = TRUE
))
summary_rows <- optimization_result$summary_rows

summary_print <- time_section("print_summary", {
  cat("Optimization objective summary\n")
  summary_print <- optimization_result$summary_print
  print(summary_print, row.names = FALSE)
  summary_print
})
if (write_summary_csv) {
  invisible(time_section("write_summary_csv", {
    ensure_parent_dir(summary_csv_file)
    write.csv(summary_print, summary_csv_file, row.names = FALSE)
    cat("summary_csv_file:", summary_csv_file, "\n")
  }))
}
if (write_portfolio_tables_csv) {
  invisible(time_section("write_portfolio_tables_csv", {
    portfolio_table_files <- write_option_portfolio_tables_csv(
      optimization_result$runs,
      portfolio_tables_csv_dir
    )
    cat("portfolio_tables_csv_dir:", portfolio_tables_csv_dir, "\n")
    cat("portfolio_table_csv_files:", length(portfolio_table_files), "\n")
  }))
}
if (write_combined_portfolio_table_csv) {
  invisible(time_section("write_combined_portfolio_table_csv", {
    write_combined_option_portfolio_table_csv(
      optimization_result$runs,
      combined_portfolio_table_csv_file
    )
    cat("combined_portfolio_table_csv_file:", combined_portfolio_table_csv_file, "\n")
  }))
}
if (write_portfolio_plot) {
  invisible(time_section("portfolio_plot", {
  ensure_parent_dir(portfolio_plot_file)
  selected_plot_run <- select_option_optimization_run(
    optimization_result$runs,
    constraint_mode = plot_constraint_mode,
    objective = plot_objective,
    utility_gamma = plot_utility_gamma
  )
  plot_option_optimization_report(
    file = portfolio_plot_file,
    selected_run = selected_plot_run,
    S0 = S0,
    r = r,
    q = q,
    T = T,
    K = K,
    type = option_type,
    market_vol = vol_quotes$mid_vol,
    expected_payoff = expected_payoff,
    rf_growth = rf_growth
  )
  cat("portfolio_plot_file:", portfolio_plot_file, "\n")
  }))
}

elapsed <- proc.time() - start_time
cat("\nSection timings\n")
section_timings_print <- section_timings
section_timings_print$elapsed_seconds <- round(section_timings_print$elapsed_seconds, 6)
print(section_timings_print, row.names = FALSE)
cat("\nelapsed_seconds:", round(elapsed["elapsed"], 6), "\n")
