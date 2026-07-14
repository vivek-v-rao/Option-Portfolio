#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Rcpp acceleration kernels for option portfolio optimization.
//
// These routines evaluate objective functions and analytic gradients used by
// xoptimize_options.r. They assume payoff matrices and price vectors have
// already been built and validated on the R side.

using namespace Rcpp;

// Net cash cost of rebalancing from initial_contracts to contracts at bid/ask.
static double trade_cost_cpp(const NumericVector& contracts,
                             const NumericVector& initial_contracts,
                             const NumericVector& bid_prices,
                             const NumericVector& ask_prices) {
  const int n = contracts.size();
  double out = 0.0;
  for (int i = 0; i < n; ++i) {
    const double trade = contracts[i] - initial_contracts[i];
    out += trade * (trade < 0.0 ? bid_prices[i] : ask_prices[i]);
  }
  return out;
}

// Piecewise derivative of trade cost with respect to one final contract count.
static double trade_cost_grad_cpp(const double contract,
                                  const double initial_contract,
                                  const double bid_price,
                                  const double ask_price) {
  return contract < initial_contract ? bid_price : ask_price;
}

// Quadratic penalty for ES/CVaR loss-cap violations across configured tails.
static double cvar_penalty_cpp(const std::vector<double>& terminal_wealth,
                               const double budget,
                               const NumericVector& tail_probs,
                               const NumericVector& max_loss,
                               const double penalty_scale) {
  const int n = static_cast<int>(terminal_wealth.size());
  const int n_tail = tail_probs.size();
  if (n == 0 || n_tail == 0) {
    return 0.0;
  }

  std::vector<double> losses(n);
  for (int i = 0; i < n; ++i) {
    losses[i] = budget - terminal_wealth[i];
  }
  std::sort(losses.begin(), losses.end(), std::greater<double>());

  double penalty_sum = 0.0;
  for (int j = 0; j < n_tail; ++j) {
    int k = static_cast<int>(std::ceil(tail_probs[j] * n));
    if (k < 1) {
      k = 1;
    } else if (k > n) {
      k = n;
    }

    double es = 0.0;
    for (int i = 0; i < k; ++i) {
      es += losses[i];
    }
    es /= static_cast<double>(k);

    const double excess = std::max(es - max_loss[j], 0.0);
    penalty_sum += excess * excess;
  }

  return penalty_scale * penalty_sum / std::max(1.0, budget * budget);
}

// Add the subgradient of the ES/CVaR penalty to an existing objective gradient.
static void add_cvar_penalty_gradient_cpp(NumericVector& grad,
                                          const NumericVector& contracts,
                                          const NumericMatrix& cvar_payoff_scenarios,
                                          const NumericVector& initial_contracts,
                                          const NumericVector& bid_prices,
                                          const NumericVector& ask_prices,
                                          const double budget,
                                          const double rf_growth,
                                          const NumericVector& tail_probs,
                                          const NumericVector& max_loss,
                                          const double penalty_scale) {
  const int n_contracts = contracts.size();
  const int n_scenarios = cvar_payoff_scenarios.nrow();
  const int n_tail = tail_probs.size();
  if (n_contracts == 0 || n_scenarios == 0 || n_tail == 0) {
    return;
  }

  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = (budget - trade_cost) * rf_growth;
  std::vector<double> trade_grad(n_contracts);
  for (int j = 0; j < n_contracts; ++j) {
    trade_grad[j] = trade_cost_grad_cpp(contracts[j], initial_contracts[j], bid_prices[j], ask_prices[j]);
  }

  std::vector<std::pair<double, int> > losses(n_scenarios);
  for (int i = 0; i < n_scenarios; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += cvar_payoff_scenarios(i, j) * contracts[j];
    }
    losses[i] = std::make_pair(budget - wealth, i);
  }
  std::sort(losses.begin(), losses.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
              return a.first > b.first;
            });

  const double denom = std::max(1.0, budget * budget);
  for (int t = 0; t < n_tail; ++t) {
    int k = static_cast<int>(std::ceil(tail_probs[t] * n_scenarios));
    if (k < 1) {
      k = 1;
    } else if (k > n_scenarios) {
      k = n_scenarios;
    }

    double es = 0.0;
    for (int i = 0; i < k; ++i) {
      es += losses[i].first;
    }
    es /= static_cast<double>(k);

    const double excess = std::max(es - max_loss[t], 0.0);
    if (excess <= 0.0) {
      continue;
    }

    const double scale = penalty_scale * 2.0 * excess / denom / static_cast<double>(k);
    for (int i = 0; i < k; ++i) {
      const int row = losses[i].second;
      for (int j = 0; j < n_contracts; ++j) {
        const double dloss = rf_growth * trade_grad[j] - cvar_payoff_scenarios(row, j);
        grad[j] -= scale * dloss;
      }
    }
  }
}

// SLSQP inequality values: linear constraints plus signed ES/CVaR excesses.
// Feasible constraints are <= 0.
// [[Rcpp::export]]
NumericVector option_constraint_values_cpp(const NumericVector& contracts,
                                           const NumericMatrix& ui,
                                           const NumericVector& ci,
                                           const NumericMatrix& cvar_payoff_scenarios,
                                           const NumericVector& initial_contracts,
                                           const NumericVector& bid_prices,
                                           const NumericVector& ask_prices,
                                           const double budget,
                                           const double rf_growth,
                                           const NumericVector& cvar_tail_probs,
                                           const NumericVector& cvar_max_loss) {
  const int n_contracts = contracts.size();
  const int n_linear = ui.nrow();
  const int n_cvar = cvar_payoff_scenarios.nrow();
  const int n_tail = cvar_tail_probs.size();
  NumericVector out(n_linear + n_tail);

  for (int i = 0; i < n_linear; ++i) {
    double value = ci[i];
    for (int j = 0; j < n_contracts; ++j) {
      value -= ui(i, j) * contracts[j];
    }
    out[i] = value;
  }

  if (n_cvar == 0 || n_tail == 0) {
    return out;
  }

  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = (budget - trade_cost) * rf_growth;
  std::vector<double> losses(n_cvar);
  for (int i = 0; i < n_cvar; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += cvar_payoff_scenarios(i, j) * contracts[j];
    }
    losses[i] = budget - wealth;
  }
  std::sort(losses.begin(), losses.end(), std::greater<double>());

  for (int t = 0; t < n_tail; ++t) {
    int k = static_cast<int>(std::ceil(cvar_tail_probs[t] * n_cvar));
    if (k < 1) {
      k = 1;
    } else if (k > n_cvar) {
      k = n_cvar;
    }

    double es = 0.0;
    for (int i = 0; i < k; ++i) {
      es += losses[i];
    }
    es /= static_cast<double>(k);
    out[n_linear + t] = es - cvar_max_loss[t];
  }

  return out;
}

// SLSQP inequality Jacobian matching option_constraint_values_cpp.
// ES/CVaR rows use the current worst-tail scenario set as a subgradient.
// [[Rcpp::export]]
NumericMatrix option_constraint_jacobian_cpp(const NumericVector& contracts,
                                             const NumericMatrix& ui,
                                             const NumericMatrix& cvar_payoff_scenarios,
                                             const NumericVector& initial_contracts,
                                             const NumericVector& bid_prices,
                                             const NumericVector& ask_prices,
                                             const double budget,
                                             const double rf_growth,
                                             const NumericVector& cvar_tail_probs) {
  const int n_contracts = contracts.size();
  const int n_linear = ui.nrow();
  const int n_cvar = cvar_payoff_scenarios.nrow();
  const int n_tail = cvar_tail_probs.size();
  NumericMatrix jac(n_linear + n_tail, n_contracts);

  for (int i = 0; i < n_linear; ++i) {
    for (int j = 0; j < n_contracts; ++j) {
      jac(i, j) = -ui(i, j);
    }
  }

  if (n_cvar == 0 || n_tail == 0) {
    return jac;
  }

  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = (budget - trade_cost) * rf_growth;
  std::vector<double> trade_grad(n_contracts);
  for (int j = 0; j < n_contracts; ++j) {
    trade_grad[j] = trade_cost_grad_cpp(contracts[j], initial_contracts[j], bid_prices[j], ask_prices[j]);
  }

  std::vector<std::pair<double, int> > losses(n_cvar);
  for (int i = 0; i < n_cvar; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += cvar_payoff_scenarios(i, j) * contracts[j];
    }
    losses[i] = std::make_pair(budget - wealth, i);
  }
  std::sort(losses.begin(), losses.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
              return a.first > b.first;
            });

  for (int t = 0; t < n_tail; ++t) {
    int k = static_cast<int>(std::ceil(cvar_tail_probs[t] * n_cvar));
    if (k < 1) {
      k = 1;
    } else if (k > n_cvar) {
      k = n_cvar;
    }

    for (int row_idx = 0; row_idx < k; ++row_idx) {
      const int row = losses[row_idx].second;
      for (int j = 0; j < n_contracts; ++j) {
        jac(n_linear + t, j) +=
          (rf_growth * trade_grad[j] - cvar_payoff_scenarios(row, j)) / static_cast<double>(k);
      }
    }
  }

  return jac;
}

// Expected CRRA/log utility of terminal wealth for unconstrained contract weights.
// Returns -Inf if any scenario has nonpositive terminal wealth.
// [[Rcpp::export]]
double expected_utility_objective_cpp(const NumericVector& contracts,
                                      const NumericMatrix& payoff_scenarios,
                                      const NumericMatrix& cvar_payoff_scenarios,
                                      const NumericVector& initial_contracts,
                                      const NumericVector& bid_prices,
                                      const NumericVector& ask_prices,
                                      const double budget,
                                      const double rf_growth,
                                      const double gamma,
                                      const NumericVector& cvar_tail_probs,
                                      const NumericVector& cvar_max_loss,
                                      const double penalty_scale = 1e6) {
  const int n_contracts = contracts.size();
  const int n_scenarios = payoff_scenarios.nrow();
  const int n_cvar = cvar_payoff_scenarios.nrow();
  const double cash = budget - trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = cash * rf_growth;

  double utility_sum = 0.0;
  for (int i = 0; i < n_scenarios; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += payoff_scenarios(i, j) * contracts[j];
    }
    if (!(wealth > 0.0) || !std::isfinite(wealth)) {
      return R_NegInf;
    }

    if (std::abs(gamma - 1.0) < 1e-10) {
      utility_sum += std::log(wealth);
    } else {
      utility_sum += std::pow(wealth, 1.0 - gamma) / (1.0 - gamma);
    }
  }

  double cvar_penalty = 0.0;
  if (n_cvar > 0 && cvar_tail_probs.size() > 0) {
    std::vector<double> cvar_wealth(n_cvar);
    for (int i = 0; i < n_cvar; ++i) {
      double wealth = cash_growth;
      for (int j = 0; j < n_contracts; ++j) {
        wealth += cvar_payoff_scenarios(i, j) * contracts[j];
      }
      cvar_wealth[i] = wealth;
    }
    cvar_penalty = cvar_penalty_cpp(cvar_wealth, budget, cvar_tail_probs, cvar_max_loss, penalty_scale);
  }

  return utility_sum / static_cast<double>(n_scenarios) - cvar_penalty;
}

// Analytic gradient of expected_utility_objective_cpp with bid/ask trade costs.
// [[Rcpp::export]]
NumericVector expected_utility_gradient_cpp(const NumericVector& contracts,
                                            const NumericMatrix& payoff_scenarios,
                                            const NumericMatrix& cvar_payoff_scenarios,
                                            const NumericVector& initial_contracts,
                                            const NumericVector& bid_prices,
                                            const NumericVector& ask_prices,
                                            const double budget,
                                            const double rf_growth,
                                            const double gamma,
                                            const NumericVector& cvar_tail_probs,
                                            const NumericVector& cvar_max_loss,
                                            const double penalty_scale = 1e6) {
  const int n_contracts = contracts.size();
  const int n_scenarios = payoff_scenarios.nrow();
  NumericVector grad(n_contracts);
  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = (budget - trade_cost) * rf_growth;
  std::vector<double> trade_grad(n_contracts);
  for (int j = 0; j < n_contracts; ++j) {
    trade_grad[j] = trade_cost_grad_cpp(contracts[j], initial_contracts[j], bid_prices[j], ask_prices[j]);
  }

  for (int i = 0; i < n_scenarios; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += payoff_scenarios(i, j) * contracts[j];
    }
    if (!(wealth > 0.0) || !std::isfinite(wealth)) {
      return NumericVector(n_contracts, NA_REAL);
    }

    const double marginal_utility = std::abs(gamma - 1.0) < 1e-10 ? 1.0 / wealth : std::pow(wealth, -gamma);
    for (int j = 0; j < n_contracts; ++j) {
      grad[j] += marginal_utility * (payoff_scenarios(i, j) - rf_growth * trade_grad[j]);
    }
  }

  for (int j = 0; j < n_contracts; ++j) {
    grad[j] /= static_cast<double>(n_scenarios);
  }
  add_cvar_penalty_gradient_cpp(
    grad, contracts, cvar_payoff_scenarios, initial_contracts, bid_prices, ask_prices,
    budget, rf_growth, cvar_tail_probs, cvar_max_loss, penalty_scale
  );
  return grad;
}

// Objective for nonnegative-terminal constrained optimization over contracts.
// Supports mean_wealth, mean_variance, and sharpe objectives.
// [[Rcpp::export]]
double nonnegative_terminal_objective_cpp(const NumericVector& contracts,
                                          const NumericVector& expected_payoff,
                                          const NumericMatrix& cov_payoff,
                                          const NumericMatrix& cvar_payoff_scenarios,
                                          const NumericVector& initial_contracts,
                                          const NumericVector& bid_prices,
                                          const NumericVector& ask_prices,
                                          const double budget,
                                          const double rf_growth,
                                          const std::string objective,
                                          const double risk_aversion,
                                          const NumericVector& cvar_tail_probs,
                                          const NumericVector& cvar_max_loss,
                                          const double tol = 1e-10,
                                          const double penalty_scale = 1e6) {
  const int n_contracts = contracts.size();
  const int n_cvar = cvar_payoff_scenarios.nrow();
  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);

  double expected_value = 0.0;
  for (int i = 0; i < n_contracts; ++i) {
    expected_value += expected_payoff[i] * contracts[i];
  }
  const double excess = expected_value - rf_growth * trade_cost;

  double variance = 0.0;
  for (int i = 0; i < n_contracts; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_contracts; ++j) {
      row_sum += cov_payoff(i, j) * contracts[j];
    }
    variance += contracts[i] * row_sum;
  }
  if (variance < 0.0 && variance > -tol) {
    variance = 0.0;
  }
  const double sd = std::sqrt(std::max(variance, 0.0));

  double cvar_penalty = 0.0;
  if (n_cvar > 0 && cvar_tail_probs.size() > 0) {
    const double cash_growth = (budget - trade_cost) * rf_growth;
    std::vector<double> cvar_wealth(n_cvar);
    for (int i = 0; i < n_cvar; ++i) {
      double wealth = cash_growth;
      for (int j = 0; j < n_contracts; ++j) {
        wealth += cvar_payoff_scenarios(i, j) * contracts[j];
      }
      cvar_wealth[i] = wealth;
    }
    cvar_penalty = cvar_penalty_cpp(cvar_wealth, budget, cvar_tail_probs, cvar_max_loss, penalty_scale);
  }

  if (objective == "mean_wealth") {
    return excess - cvar_penalty;
  }
  if (objective == "mean_variance") {
    return excess / budget - risk_aversion * variance / (budget * budget) - cvar_penalty;
  }
  if (!(std::isfinite(sd)) || sd <= tol) {
    return R_NegInf;
  }
  if (objective == "sharpe") {
    return excess / sd - cvar_penalty;
  }

  return R_NegInf;
}

// Analytic gradient of nonnegative_terminal_objective_cpp.
// [[Rcpp::export]]
NumericVector nonnegative_terminal_gradient_cpp(const NumericVector& contracts,
                                                const NumericVector& expected_payoff,
                                                const NumericMatrix& cov_payoff,
                                                const NumericMatrix& cvar_payoff_scenarios,
                                                const NumericVector& initial_contracts,
                                                const NumericVector& bid_prices,
                                                const NumericVector& ask_prices,
                                                const double budget,
                                                const double rf_growth,
                                                const std::string objective,
                                                const double risk_aversion,
                                                const NumericVector& cvar_tail_probs,
                                                const NumericVector& cvar_max_loss,
                                                const double tol = 1e-10,
                                                const double penalty_scale = 1e6) {
  const int n_contracts = contracts.size();
  NumericVector grad(n_contracts);
  std::vector<double> trade_grad(n_contracts);
  for (int i = 0; i < n_contracts; ++i) {
    trade_grad[i] = trade_cost_grad_cpp(contracts[i], initial_contracts[i], bid_prices[i], ask_prices[i]);
  }

  double expected_value = 0.0;
  for (int i = 0; i < n_contracts; ++i) {
    expected_value += expected_payoff[i] * contracts[i];
  }
  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double excess = expected_value - rf_growth * trade_cost;

  NumericVector cov_contracts(n_contracts);
  double variance = 0.0;
  for (int i = 0; i < n_contracts; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_contracts; ++j) {
      row_sum += cov_payoff(i, j) * contracts[j];
    }
    cov_contracts[i] = row_sum;
    variance += contracts[i] * row_sum;
  }
  if (variance < 0.0 && variance > -tol) {
    variance = 0.0;
  }
  const double sd = std::sqrt(std::max(variance, 0.0));

  for (int i = 0; i < n_contracts; ++i) {
    const double dexcess = expected_payoff[i] - rf_growth * trade_grad[i];
    if (objective == "mean_wealth") {
      grad[i] = dexcess;
    } else if (objective == "mean_variance") {
      grad[i] = dexcess / budget - 2.0 * risk_aversion * cov_contracts[i] / (budget * budget);
    } else if (objective == "sharpe") {
      if (!(std::isfinite(sd)) || sd <= tol) {
        return NumericVector(n_contracts, NA_REAL);
      }
      grad[i] = dexcess / sd - excess * cov_contracts[i] / (sd * sd * sd);
    } else {
      return NumericVector(n_contracts, NA_REAL);
    }
  }

  add_cvar_penalty_gradient_cpp(
    grad, contracts, cvar_payoff_scenarios, initial_contracts, bid_prices, ask_prices,
    budget, rf_growth, cvar_tail_probs, cvar_max_loss, penalty_scale
  );
  return grad;
}

// Expected utility for long-only portfolios parameterized by softmax weights.
// The optional cash weight is included when force_full_investment is false.
// [[Rcpp::export]]
double long_only_expected_utility_objective_cpp(const NumericVector& theta,
                                                const NumericMatrix& payoff_scenarios,
                                                const NumericVector& prices,
                                                const NumericVector& initial_contracts,
                                                const NumericVector& bid_prices,
                                                const NumericVector& ask_prices,
                                                const double budget,
                                                const double rf_growth,
                                                const double gamma,
                                                const double max_weight,
                                                const bool force_full_investment) {
  const int n_contracts = prices.size();
  const int n_theta = theta.size();
  const int n_scenarios = payoff_scenarios.nrow();
  const int softmax_len = force_full_investment ? n_contracts : n_contracts + 1;
  if (n_theta != softmax_len) {
    return R_NegInf;
  }

  double theta_max = theta[0];
  for (int i = 1; i < n_theta; ++i) {
    if (theta[i] > theta_max) {
      theta_max = theta[i];
    }
  }

  std::vector<double> exp_theta(n_theta);
  double exp_sum = 0.0;
  for (int i = 0; i < n_theta; ++i) {
    exp_theta[i] = std::exp(theta[i] - theta_max);
    exp_sum += exp_theta[i];
  }
  if (!(exp_sum > 0.0) || !std::isfinite(exp_sum)) {
    return R_NegInf;
  }

  NumericVector contracts(n_contracts);
  for (int i = 0; i < n_contracts; ++i) {
    const double weight = max_weight * exp_theta[i] / exp_sum;
    contracts[i] = budget * weight / prices[i];
  }

  const double cash = budget - trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = cash * rf_growth;

  double utility_sum = 0.0;
  for (int i = 0; i < n_scenarios; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += payoff_scenarios(i, j) * contracts[j];
    }
    if (!(wealth > 0.0) || !std::isfinite(wealth)) {
      return R_NegInf;
    }

    if (std::abs(gamma - 1.0) < 1e-10) {
      utility_sum += std::log(wealth);
    } else {
      utility_sum += std::pow(wealth, 1.0 - gamma) / (1.0 - gamma);
    }
  }

  return utility_sum / static_cast<double>(n_scenarios);
}

// Gradient of long_only_expected_utility_objective_cpp with respect to theta.
// [[Rcpp::export]]
NumericVector long_only_expected_utility_gradient_cpp(const NumericVector& theta,
                                                      const NumericMatrix& payoff_scenarios,
                                                      const NumericVector& prices,
                                                      const NumericVector& initial_contracts,
                                                      const NumericVector& bid_prices,
                                                      const NumericVector& ask_prices,
                                                      const double budget,
                                                      const double rf_growth,
                                                      const double gamma,
                                                      const double max_weight,
                                                      const bool force_full_investment) {
  const int n_contracts = prices.size();
  const int n_theta = theta.size();
  const int n_scenarios = payoff_scenarios.nrow();
  const int softmax_len = force_full_investment ? n_contracts : n_contracts + 1;
  if (n_theta != softmax_len) {
    return NumericVector(n_theta, NA_REAL);
  }

  double theta_max = theta[0];
  for (int i = 1; i < n_theta; ++i) {
    if (theta[i] > theta_max) {
      theta_max = theta[i];
    }
  }

  std::vector<double> prob(n_theta);
  double exp_sum = 0.0;
  for (int i = 0; i < n_theta; ++i) {
    prob[i] = std::exp(theta[i] - theta_max);
    exp_sum += prob[i];
  }
  if (!(exp_sum > 0.0) || !std::isfinite(exp_sum)) {
    return NumericVector(n_theta, NA_REAL);
  }
  for (int i = 0; i < n_theta; ++i) {
    prob[i] /= exp_sum;
  }

  NumericVector contracts(n_contracts);
  for (int j = 0; j < n_contracts; ++j) {
    contracts[j] = budget * max_weight * prob[j] / prices[j];
  }

  NumericVector contract_grad(n_contracts);
  const double trade_cost = trade_cost_cpp(contracts, initial_contracts, bid_prices, ask_prices);
  const double cash_growth = (budget - trade_cost) * rf_growth;
  std::vector<double> trade_grad(n_contracts);
  for (int j = 0; j < n_contracts; ++j) {
    trade_grad[j] = trade_cost_grad_cpp(contracts[j], initial_contracts[j], bid_prices[j], ask_prices[j]);
  }

  for (int i = 0; i < n_scenarios; ++i) {
    double wealth = cash_growth;
    for (int j = 0; j < n_contracts; ++j) {
      wealth += payoff_scenarios(i, j) * contracts[j];
    }
    if (!(wealth > 0.0) || !std::isfinite(wealth)) {
      return NumericVector(n_theta, NA_REAL);
    }

    const double marginal_utility = std::abs(gamma - 1.0) < 1e-10 ? 1.0 / wealth : std::pow(wealth, -gamma);
    for (int j = 0; j < n_contracts; ++j) {
      contract_grad[j] += marginal_utility * (payoff_scenarios(i, j) - rf_growth * trade_grad[j]);
    }
  }
  for (int j = 0; j < n_contracts; ++j) {
    contract_grad[j] /= static_cast<double>(n_scenarios);
  }

  double weighted_contract_grad = 0.0;
  for (int j = 0; j < n_contracts; ++j) {
    weighted_contract_grad += contract_grad[j] * contracts[j];
  }

  NumericVector theta_grad(n_theta);
  for (int k = 0; k < n_theta; ++k) {
    theta_grad[k] = -prob[k] * weighted_contract_grad;
    if (k < n_contracts) {
      theta_grad[k] += contract_grad[k] * contracts[k];
    }
  }
  return theta_grad;
}
