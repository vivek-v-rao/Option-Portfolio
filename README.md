# Optimal Options

R code for optimizing option portfolios under multiple objectives and constraints. The main driver is single-underlying and single-expiration; prototype drivers cover named structures, multiple tenors, and multiple stocks.

## Files

- `xoptimize_options.r`: main example driver. Edit top-level parameters here.
- `xoptimize_option_structures.r`: named-structure search driver for straddles, butterflies, iron condors, and related templates.
- `xoptimize_options_multitenor.r`: prototype multi-expiration optimizer that marks all instruments at the first expiry.
- `xoptimize_options_multistock.r`: prototype single-expiration, multi-underlying optimizer using multivariate lognormal terminal prices.
- `xoptimize_options_quick.r`: fast smoke test that runs the main driver with a smaller strike grid and fewer CVaR scenarios.
- `option_stats.r`: reusable pricing, payoff moment, transaction cost, constraint, optimizer, reporting, and utility functions.
- `configs/`: JSON run configurations intended to be shared by the R and future Python drivers.
- `outputs/`: generated CSV and plot output directory; ignored by Git.

## Model

The main `xoptimize_options.r` driver assumes one underlying stock and one option expiration. Prototype drivers extend the workflow to named structures, multiple tenors, and multiple stocks. In the main driver, the terminal stock price model is selected with `terminal_model`:

- `lognormal`: default GBM-style terminal stock price.
- `logistic`: logistic distribution for log returns, with scale chosen to match `realized_sigma`.
- `hyperbolic_secant`: hyperbolic secant distribution for log returns, with scale chosen to match `realized_sigma`.
- `symmetric_hyperbolic`: symmetric hyperbolic distribution for log returns, with `hyperbolic_alpha` controlling tail thickness.
- `hyperbolic`: asymmetric hyperbolic distribution for log returns, with `hyperbolic_alpha` controlling tail thickness and `hyperbolic_beta` controlling skew.
- `generalized_hyperbolic`: generalized hyperbolic distribution for log returns, with `gh_lambda`, `gh_alpha`, and `gh_beta`.
- `normal_inverse_gaussian`: NIG distribution for log returns, equivalent to generalized hyperbolic with lambda -1/2.
- `variance_gamma`: variance-gamma distribution for log returns, with `vg_shape` and `vg_beta`.
- `lognormal_mixture`: mixture of GBM-style lognormal terminal regimes.
- `normal`: Bachelier-style price-space terminal distribution.

For price-space distributions, set `terminal_floor` to control limited liability:

- `terminal_floor <- 0.0`: values below zero are collapsed into a point mass at zero.
- `terminal_floor <- NA_real_`: no floor; negative terminal prices are allowed.

The default lognormal setup uses user-specified drift and realized volatility:

```r
S0 <- 100
mu <- 0.03
realized_sigma <- 0.20
T <- 1.0
r <- 0.03
q <- 0.00
terminal_model <- "lognormal"
terminal_floor <- NA_real_
```

For `lognormal_mixture`, set vectors of equal length:

```r
mixture_weights <- c(0.85, 0.15)
mixture_mu <- c(0.06, -0.25)
mixture_sigma <- c(0.16, 0.45)
```

Option prices are still generated from the configured `implied_sigma`, so a mixture scenario can compare flat-vol market prices against a richer assumed true terminal distribution.

For `logistic`, the code matches `E[S_T] = S0 * exp(mu*T)` and `Var(log(S_T)) = realized_sigma^2*T`. It requires the logistic log-return scale to be below 1, which ensures finite expected stock and call payoffs.

For `hyperbolic_secant`, the code uses the standard hyperbolic secant distribution with variance 1, then scales it so `Var(log(S_T)) = realized_sigma^2*T` and adjusts location so `E[S_T] = S0 * exp(mu*T)`. It requires `realized_sigma * sqrt(T) < pi/2` for finite expected stock and call payoffs.

For `symmetric_hyperbolic`, the code uses a symmetric hyperbolic log-return distribution. `hyperbolic_alpha` controls exponential tail decay and must be greater than 1 for finite call prices. The scale parameter is solved so `Var(log(S_T)) = realized_sigma^2*T`, and `hyperbolic_alpha` must be large enough for that variance match.

For `hyperbolic`, `hyperbolic_beta` adds asymmetry. Negative values create heavier downside tails in log returns. The code requires `hyperbolic_alpha > abs(hyperbolic_beta + 1)` so expected stock and call payoffs are finite, and solves the scale parameter to match `Var(log(S_T)) = realized_sigma^2*T`.

For `generalized_hyperbolic`, `normal_inverse_gaussian`, and `variance_gamma`, terminal quantiles are computed numerically from the configured density. The location parameter is adjusted to match `E[S_T] = S0 * exp(mu*T)`, and the scale/tail parameter is solved to match `Var(log(S_T)) = realized_sigma^2*T` when needed.

The optimizer is a static terminal-payoff model. It assumes the option, stock, and cash positions are established at the initial bid/ask prices and then held to expiration. It does not simulate delta hedging, interim rebalancing, path-dependent P&L, hedge transaction costs, or implied-volatility marks before expiration. A delta-hedged vol-arbitrage model would require a full return path process, hedge schedule, and hedge cost model rather than only a terminal distribution.

Separately, the program can report short-horizon mark-to-market VaR/ES by repricing the portfolio after a one-period return shock. This diagnostic currently supports `var_return_model <- "normal_log"` and `var_return_model <- "normal_simple"`, uses sticky-strike mid vols, and is not part of the continuous optimizer constraint set. The multi-stock prototype uses correlated one-period return shocks from the configured covariance matrix.

It can also report simulated delta-hedged terminal wealth diagnostics. Set:

```r
delta_hedge_steps <- 0L
delta_hedge_paths <- 1000L
delta_hedge_mu <- NA_real_
delta_hedge_sigma <- NA_real_
delta_hedge_seed <- NA_integer_
delta_hedge_stock_transaction_cost <- 0.0
```

`delta_hedge_steps <- 0L` is the default and means no delta hedging. Positive values simulate that many equally spaced hedge dates before expiration under GBM with sticky-strike implied vols for hedge deltas. `delta_hedge_mu` and `delta_hedge_sigma` default to `mu` and `realized_sigma`. The hedged simulation is currently a reporting layer, not an optimizer objective or native constraint.

The instrument universe contains:

- the stock, shown as `type == "call"` and `strike == 0`
- calls at `base_strikes`
- puts at the same `base_strikes`

Calls and puts use the same strike-level `implied_sigma` vector. Bid/ask transaction costs can be modeled with bid/ask implied volatilities, where the vol spread widens quadratically in log-moneyness. The stock row has its own bid/ask spread.

## Objectives

Set `optimization_objectives` in `xoptimize_options.r`.

Currently supported:

- `sharpe`
- `mean_variance`
- `adjusted_sharpe` or `sharpe_adj`
- `dependency_penalty`
- `expected_utility`

For `expected_utility`, `risk_aversion_utility` may be a vector, and the program runs one optimization per value.

`dependency_penalty` estimates a conditional-payout dependency matrix from the terminal payoff scenarios, repairs it to be positive semidefinite by flooring eigenvalues, and maximizes edge minus `risk_aversion * w' Lambda w`. This follows the idea in `portfolio_optimisation_with_options.pdf`: options whose positive payoffs tend to occur in the same scenarios are penalized as dependent.

## Constraints

Set `constraint_modes` in `xoptimize_options.r`.

Currently supported:

- `long_only`
- `long_short`
- `nonnegative_terminal`

Use `long_short` when the intended trade is to buy options with cheap implied vol and sell options with rich implied vol. It optimizes signed instrument weights by splitting each instrument into long and short sleeves; ask prices are used for long trades and bid prices for short trades. `max_invested_weight` is interpreted as a gross signed exposure cap, and the summary reports both net `invested_weight` and `gross_position_weight`.

The `nonnegative_terminal` mode allows short option structures, but scales the risky overlay toward cash to satisfy terminal-wealth constraints. It can also enforce:

- `min_terminal_wealth`
- multiple CVaR/ES loss caps via `cvar_constraints`
- maximum absolute position weight via `max_invested_weight`
- optional portfolio delta bounds via `portfolio_delta_bounds`

Use `portfolio_delta_bounds <- c(lower, upper)` to limit total portfolio delta in share units, computed from mid-vol Black-Scholes deltas. For example:

```r
portfolio_delta_bounds <- c(-0.05, 0.05)
```

The default `c(NA_real_, NA_real_)` leaves delta unconstrained. Delta bounds are native linear constraints for `nonnegative_terminal`. For `long_only`, the optimizer scales the risky overlay toward the initial portfolio if the unconstrained solution breaches the band; tight bounds may therefore leave more cash invested than an unconstrained long-only run.

## Transaction Costs And Initial Positions

Set `use_bid_ask_vols <- TRUE` to trade options at bid/ask prices generated from the quadratic vol-spread model.

Set `initial_contracts` to represent an existing portfolio. Transaction costs are charged only on trades from `initial_contracts` to optimized final contracts.

## Integer Contracts

By default the optimizer reports continuous contract quantities. Set:

```r
integer_contracts <- TRUE
integer_rounding_neighborhood <- 1L
integer_max_search_instruments <- 6L
```

to round trade contracts to integers and search nearby integer portfolios. The continuous solution is still used as the starting point. The final summary reports `continuous_obj_value`, signed `integer_obj_delta`, nonnegative `integer_obj_loss`, `integer_candidate_count`, `integer_feasible_count`, and `integer_used_initial_fallback` so the cost and quality of the integer repair are visible.

## Simple Portfolios

To cap the number of option instruments with active positions, set:

```r
max_active_option_positions <- 3L
prune_positions_by <- "position_weight" # or "trade_weight", "contracts"
prune_repair_constraints <- TRUE
```

This is a post-optimization pruning heuristic, not a true cardinality optimizer. It keeps the largest option positions, excludes the stock row from the count, sets smaller option positions to zero, optionally repairs constraints by scaling toward the initial portfolio, and then runs integer repair if enabled. The final summary reports `pruned_obj_value`, `prune_obj_delta`, `prune_obj_loss`, `active_option_positions`, and `pruned_option_positions`.

## Named Structures

Use `xoptimize_option_structures.r` when the trader wants a specific payoff template and only wants the program to choose strikes:

```powershell
Rscript .\xoptimize_option_structures.r
```

The driver generates candidates from `base_strikes`, evaluates each fixed structure with the same pricing, payoff, Greek, MTM VaR/ES, and optional delta-hedged diagnostics, then ranks candidates by `structure_objective`.

Supported structure templates include:

- `straddle`: long call and long put at one strike.
- `strangle`: long put at lower strike and long call at higher strike.
- `call_spread`: long lower-strike call and short higher-strike call.
- `put_spread`: short lower-strike put and long higher-strike put.
- `butterfly`: long/short/long call butterfly.
- `iron_condor`: long put wing, short put, short call, long call wing.

Set:

```r
structure_types <- c("straddle", "butterfly", "iron_condor")
structure_units <- 1.0
structure_objective <- "sharpe"
structure_top_n <- 10L
```

This is a strike/template search, not a free portfolio optimizer.

## Multiple Tenors

Use `xoptimize_options_multitenor.r` for the first multi-expiration workflow:

```powershell
Rscript .\xoptimize_options_multitenor.r
```

The prototype optimizes to `horizon <- min(expiries)`. Options expiring at the horizon settle to intrinsic payoff; later-tenor options are marked at the horizon with Black-Scholes using remaining maturity and sticky-strike implied vols. The reusable helpers are:

- `multitenor_option_value_matrix()`
- `multitenor_price_vec()`
- `multitenor_greeks_table()`
- `multitenor_horizon_inputs()`

The default driver is intentionally small and long-only so it runs quickly. Broader strike grids, more expiries, nonnegative-terminal constraints, and alternative forward-vol rules should be added incrementally. For multi-tenor books, vega should normally be reported by expiry bucket. A single total vega is meaningful only for a specified parallel implied-volatility shift across tenors; calendar trades should keep signed vega by tenor visible.

## Multiple Stocks

Use `xoptimize_options_multistock.r` for the first single-expiration, multi-underlying workflow:

```powershell
Rscript .\xoptimize_options_multistock.r
```

The prototype simulates terminal prices from a multivariate lognormal distribution with user-specified drifts, volatilities, and correlation matrix. It builds one combined scenario payoff matrix for all stock and option instruments, then reuses the existing optimizer and reporting code.

The first version is intentionally narrow: it compares long-only and long/short constraints, reports aggregate greeks, and reports one-period correlated mark-to-market VaR/ES. For multi-stock books, delta, gamma, and vega should normally be reported by underlying; dollar delta can be summed as a rough directional exposure, but raw vega across different stocks should not be netted unless a common volatility factor model is specified. Natural next extensions are per-underlying delta/gamma/vega summaries, per-underlying greek bounds, and scenario-based nonnegative wealth constraints across the joint terminal grid.

## Running

From `C:\rcode\optimal_options`:

```powershell
Rscript .\xoptimize_options_quick.r
```

Use the quick script after refactors. It runs both constraint modes and the main objective types with a smaller setup.

For the full example:

```powershell
Rscript .\xoptimize_options.r
```

The full run can be much slower, especially for expected utility under the nonnegative terminal constraint.

You can also run from a JSON config file:

```powershell
Rscript .\xoptimize_options.r --config .\configs\quick.json
```

JSON config loading requires the R package `jsonlite`.

Config files have this shape:

```json
{
  "scenario": "quick",
  "config": {
    "optimization_objectives": ["mean_variance", "dependency_penalty"],
    "constraint_modes": ["long_only"],
    "write_summary_csv": true,
    "summary_csv_file": "outputs/example_summary.csv"
  }
}
```

The `scenario` field selects a preset. Values under `config` override that preset. Command-line `--scenario name` overrides the config-file scenario. Existing source-time variables still work: `option_scenario`, `option_config_file`, and `option_config_overrides`.

## Scenario Presets

Set `option_scenario` before sourcing `xoptimize_options.r`, or edit the default `scenario <- "full"` near the top of the file.

Available presets:

- `full`: default full example.
- `quick`: smaller strike grid and fewer CVaR scenarios for smoke tests.
- `no_costs`: disables option and stock bid/ask transaction costs.
- `long_only`: runs only the long-only constraint mode.
- `long_short`: runs signed long/short mean-variance and Sharpe optimizations.
- `logistic`: uses a logistic log-return terminal distribution.
- `hyperbolic_secant`: uses a hyperbolic secant log-return terminal distribution.
- `symmetric_hyperbolic`: uses a symmetric hyperbolic log-return terminal distribution.
- `hyperbolic`: uses an asymmetric hyperbolic log-return terminal distribution.
- `generalized_hyperbolic`: uses a generalized hyperbolic log-return terminal distribution.
- `normal_inverse_gaussian`: uses an NIG log-return terminal distribution.
- `variance_gamma`: uses a variance-gamma log-return terminal distribution.
- `normal_floor`: uses the normal price-space terminal model with `terminal_floor <- 0.0`.
- `normal_unbounded`: uses the normal price-space terminal model with `terminal_floor <- NA_real_`.
- `crash_mixture`: uses a two-regime lognormal mixture with flat implied vol option prices.
- `tail_risk`: focuses on nonnegative terminal wealth with tighter minimum wealth and ES constraints.

Preset values are applied after the base `config` list. Explicit `option_config_overrides` are applied last, so ad hoc experiments can still override preset values.

## Output

For each objective and constraint mode, the program prints an option table with:

- type, strike, bid/mid/ask vols and prices
- trade price
- payoff per dollar, edge, Sharpe, adjusted Sharpe
- Black-Scholes delta, gamma, and vega computed at mid vol
- correlation with the optimized portfolio
- initial contracts, trade contracts, trade cost, trade weight
- position value, position weight, final contracts

It also prints portfolio statistics:

- mean wealth
- standard deviation of wealth
- Sharpe and adjusted Sharpe
- skew
- excess kurtosis
- minimum terminal wealth
- tail slope
- portfolio delta, gamma, vega, delta dollars, gamma dollars for a 1% spot move, and vega for a 1 vol point move
- delta bound violation
- ES loss columns for configured CVaR tail probabilities
- short-horizon mark-to-market wealth, P&L mean/sd, VaR, and ES when `report_mtm_var` is enabled
- simulated delta-hedged wealth, P&L mean/sd, Sharpe, skew, and excess kurtosis when `delta_hedge_steps > 0`

The stock row reports delta `exp(-q*T)`, gamma `0`, and vega `0`. Option vega is per 1.00 volatility unit; `vega_1pct` is scaled to one volatility point. For single-underlying, single-expiry books, scalar portfolio delta/gamma/vega are natural summaries. For multi-stock or multi-tenor books, treat scalar Greeks as rough diagnostics and prefer bucketed reports by `underlying`, `expiry`, or `(underlying, expiry)`.

Set:

```r
report_mtm_var <- TRUE
var_horizon_days <- 1.0
var_trading_days <- 252.0
var_n_scenarios <- 1001
var_return_model <- "normal_log" # or "normal_simple"
var_conf_levels <- c(0.95, 0.99)
```

`var_return_mu` and `var_return_sigma` default to `mu` and `realized_sigma`. The reported `var_95`, `es_95`, `var_99`, and `es_99` columns are positive P&L losses over the configured horizon.

Optional interim MTM loss caps can be reported with:

```r
mtm_var_constraints <- data.frame(conf_level = c(0.95), max_loss = c(50.0))
mtm_es_constraints <- data.frame(conf_level = c(0.99), max_loss = c(100.0))
```

The summary then includes `max_mtm_var_violation`, `max_mtm_es_violation`, and per-level violation columns such as `mtm_var_violation_95`. These constraints are enforced during integer contract repair. Continuous optimization currently reports the violations rather than treating them as native optimizer constraints.

The final `Optimization objective summary` table compares all objective/constraint runs. Set:

```r
write_summary_csv <- TRUE
summary_csv_file <- "optimization_objective_summary.csv"
```

to write that summary table to CSV.

When `dependency_penalty` is enabled, the summary includes the dependency objective, penalty value, raw/repaired dependency-matrix eigenvalue diagnostics, condition number, and whether PSD repair was applied.

Set:

```r
write_portfolio_tables_csv <- TRUE
portfolio_tables_csv_dir <- "portfolio_tables"
```

to write one raw, unrounded option table CSV per objective/constraint run. Set:

```r
write_combined_portfolio_table_csv <- TRUE
combined_portfolio_table_csv_file <- "portfolio_tables/portfolio_tables_combined.csv"
```

to write one long-format CSV stacking all option-table rows with `run_id`, `constraint_mode`, `objective`, and `utility_gamma`.

Set:

```r
write_portfolio_plot <- TRUE
portfolio_plot_file <- "option_portfolio_report.png"
plot_constraint_mode <- "long_only"
plot_objective <- "mean_variance"
plot_utility_gamma <- NA_real_
```

to write a three-panel PNG report for one selected run. The report shows market versus model-implied true volatility, payoff edge by strike, and recommended position weights. For example:

```powershell
Rscript -e "option_scenario <- 'crash_mixture'; option_config_overrides <- list(write_portfolio_plot=TRUE, portfolio_plot_file='crash_mixture_plot.png', plot_constraint_mode='long_only', plot_objective='mean_variance'); source('xoptimize_options.r')"
```

## Development Notes

Keep reusable functions in `option_stats.r`; keep `xoptimize_options.r` focused on parameter setup and printing.

Scenario inputs are collected in the `config` list near the top of `xoptimize_options.r`. Alternate drivers can set `option_scenario` and/or `option_config_overrides` before sourcing `xoptimize_options.r`; `xoptimize_options_quick.r` uses this path for smoke tests.

Optional C++ kernels in `option_kernels.cpp` accelerate expected utility, constrained objectives, gradients, and SLSQP constraint callbacks when enabled with:

```r
options(option_stats_use_cpp_kernels = TRUE)
```

For constrained SLSQP runs, `constrained_optimizer_max_starts` can reduce repeated starts for Sharpe and mean-variance objectives. The nonnegative-terminal `dependency_penalty` objective currently uses the base-R constrained optimizer path. Expected utility uses the separate `nonnegative_expected_utility_max_starts` setting.

Before committing changes, run:

```powershell
Rscript .\tests\test_config_quick.r
Rscript .\xoptimize_options_quick.r
```

The config regression test runs the JSON driver path, writes test CSVs under `outputs/tests/`, and checks summary/portfolio schemas plus basic feasibility invariants.

For changes touching optimizer dispatch, payoff evaluation, or constraints, also run:

```powershell
Rscript .\xoptimize_options.r
```
