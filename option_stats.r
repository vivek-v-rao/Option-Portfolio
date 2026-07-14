# Standard normal CDF
Phi <- function(x) pnorm(x)

option_kernels_loaded <- local({
  loaded <- FALSE
  kernel_file <- file.path(getwd(), "option_kernels.cpp")
  use_cpp <- isTRUE(getOption("option_stats_use_cpp_kernels", FALSE))
  if (use_cpp && file.exists(kernel_file) && requireNamespace("Rcpp", quietly = TRUE)) {
    loaded <- tryCatch({
      old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
      makevars_file <- tempfile("option_kernels_Makevars")
      writeLines(
        c(
          "CXX = g++",
          "CXX11 = g++",
          "CXX11STD = -std=gnu++11",
          "CXX14 = g++",
          "CXX14STD = -std=gnu++14",
          "CXX17 = g++",
          "CXX17STD = -std=gnu++17",
          "SHLIB_CXXLD = g++",
          "SHLIB_CXX11LD = g++",
          "SHLIB_CXX14LD = g++",
          "SHLIB_CXX17LD = g++"
        ),
        makevars_file
      )
      Sys.setenv(R_MAKEVARS_USER = makevars_file)
      on.exit({
        if (is.na(old_makevars)) {
          Sys.unsetenv("R_MAKEVARS_USER")
        } else {
          Sys.setenv(R_MAKEVARS_USER = old_makevars)
        }
        if (file.exists(makevars_file)) {
          unlink(makevars_file)
        }
      }, add = TRUE)
      Rcpp::sourceCpp(kernel_file, rebuild = FALSE, verbose = FALSE)
      exists("expected_utility_objective_cpp", mode = "function", inherits = TRUE)
    }, error = function(e) {
      warning("Could not load option_kernels.cpp; using R fallback: ", conditionMessage(e))
      FALSE
    })
  }
  loaded
})

# d_n for lognormal truncated moments
ln_d <- function(m, v, K, n) {
  (m + n * v - log(K)) / sqrt(v)
}

# E[S_T^n * 1{S_T > K}], where log(S_T) ~ N(m, v)
ln_trunc_moment <- function(m, v, K, n) {
  exp(n * m + 0.5 * n^2 * v) * Phi(ln_d(m, v, K, n))
}

# E[S_T^n * 1{lower < S_T <= upper}]
ln_interval_moment <- function(m, v, lower, upper, n) {
  if (upper == Inf) {
    upper_moment <- 0.0
  } else {
    upper_moment <- ln_trunc_moment(m, v, upper, n)
  }

  if (lower <= 0) {
    lower_moment <- exp(n * m + 0.5 * n^2 * v)
  } else {
    lower_moment <- ln_trunc_moment(m, v, lower, n)
  }

  lower_moment - upper_moment
}

# E[((S_T - K)^+)^n]
call_payoff_raw_moment <- function(m, v, K, n) {
  out <- 0.0

  for (j in 0:n) {
    out <- out +
      choose(n, j) *
      (-K)^(n - j) *
      ln_trunc_moment(m, v, K, j)
  }

  out
}

# E[(S_T - K)^+]
call_payoff_mean <- function(m, v, K) {
  call_payoff_raw_moment(m, v, K, 1)
}

# Var[(S_T - K)^+]
call_payoff_variance <- function(m, v, K) {
  r1 <- call_payoff_raw_moment(m, v, K, 1)
  r2 <- call_payoff_raw_moment(m, v, K, 2)

  r2 - r1^2
}

# Skewness of (S_T - K)^+
call_payoff_skew <- function(m, v, K) {
  r1 <- call_payoff_raw_moment(m, v, K, 1)
  r2 <- call_payoff_raw_moment(m, v, K, 2)
  r3 <- call_payoff_raw_moment(m, v, K, 3)

  mu2 <- r2 - r1^2
  mu3 <- r3 - 3 * r1 * r2 + 2 * r1^3

  mu3 / mu2^(3 / 2)
}

# Excess kurtosis of (S_T - K)^+
call_payoff_excess_kurtosis <- function(m, v, K) {
  r1 <- call_payoff_raw_moment(m, v, K, 1)
  r2 <- call_payoff_raw_moment(m, v, K, 2)
  r3 <- call_payoff_raw_moment(m, v, K, 3)
  r4 <- call_payoff_raw_moment(m, v, K, 4)

  mu2 <- r2 - r1^2
  mu4 <- r4 - 4 * r1 * r3 + 6 * r1^2 * r2 - 3 * r1^4

  mu4 / mu2^2 - 3
}

# Convert S0, mu, sigma, T to lognormal parameters m, v
ln_params_from_gbm <- function(S0, mu, sigma, T) {
  list(
    m = log(S0) + (mu - 0.5 * sigma^2) * T,
    v = sigma^2 * T
  )
}

# E[|X - E[X]|], where log(X) ~ N(m, v)
lognormal_mean_abs_deviation <- function(m, v) {
  sigma <- sqrt(v)
  mean_x <- exp(m + 0.5 * v)
  2.0 * mean_x * (2.0 * Phi(0.5 * sigma) - 1.0)
}

lognormal_mean_average_deviation <- lognormal_mean_abs_deviation

apply_terminal_floor <- function(S, terminal_floor = NA_real_) {
  if (!is.na(terminal_floor)) {
    return(pmax(S, terminal_floor))
  }
  S
}

normalize_mixture_weights <- function(weights) {
  if (any(!is.finite(weights)) || any(weights < 0.0) || sum(weights) <= 0.0) {
    stop("mixture_weights must be nonnegative finite values with positive sum")
  }
  weights / sum(weights)
}

lognormal_mixture_params <- function(S0, mu, sigma, T) {
  stopifnot(length(mu) == length(sigma))
  list(
    m = log(S0) + (mu - 0.5 * sigma^2) * T,
    v = sigma^2 * T
  )
}

logistic_logreturn_params <- function(S0, mu, sigma, T) {
  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, scale = 0.0, v = v))
  }

  scale <- sqrt(3.0 * v) / pi
  if (scale >= 1.0) {
    stop("logistic terminal_model requires sqrt(3 * sigma^2 * T) / pi < 1 for finite call prices")
  }

  log_mgf_centered_at_one <- log(pi * scale) - log(sin(pi * scale))
  list(
    location = log(S0) + mu * T - log_mgf_centered_at_one,
    scale = scale,
    v = v
  )
}

hyperbolic_secant_logreturn_params <- function(S0, mu, sigma, T) {
  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, scale = 0.0, v = v))
  }

  scale <- sqrt(v)
  if (scale >= 0.5 * pi) {
    stop("hyperbolic_secant terminal_model requires sigma * sqrt(T) < pi / 2 for finite call prices")
  }

  list(
    location = log(S0) + mu * T + log(cos(scale)),
    scale = scale,
    v = v
  )
}

qhyperbolic_secant <- function(p, location = 0.0, scale = 1.0) {
  location + scale * (2.0 / pi) * log(tan(0.5 * pi * p))
}

symmetric_hyperbolic_variance_from_z <- function(z, alpha) {
  z / alpha^2 * besselK(z, 2) / besselK(z, 1)
}

symmetric_hyperbolic_solve_delta <- function(alpha, v) {
  min_v <- 2.0 / alpha^2
  if (v <= min_v) {
    stop("symmetric_hyperbolic terminal_model requires sigma^2 * T > 2 / hyperbolic_alpha^2; increase hyperbolic_alpha")
  }

  f <- function(z) symmetric_hyperbolic_variance_from_z(z, alpha) - v
  lo <- 1e-6
  hi <- max(1.0, alpha^2 * v)
  while (f(hi) < 0.0) {
    hi <- 2.0 * hi
  }
  uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root / alpha
}

symmetric_hyperbolic_logreturn_params <- function(S0, mu, sigma, T, alpha = 10.0) {
  if (!is.finite(alpha) || alpha <= 1.0) {
    stop("symmetric_hyperbolic terminal_model requires hyperbolic_alpha > 1 for finite call prices")
  }

  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, alpha = alpha, delta = 0.0, v = v))
  }

  delta <- symmetric_hyperbolic_solve_delta(alpha, v)
  log_mgf_centered_at_one <-
    log(alpha) - 0.5 * log(alpha^2 - 1.0) +
    log(besselK(delta * sqrt(alpha^2 - 1.0), 1)) -
    log(besselK(delta * alpha, 1))
  list(
    location = log(S0) + mu * T - log_mgf_centered_at_one,
    alpha = alpha,
    delta = delta,
    v = v
  )
}

qsymmetric_hyperbolic <- function(p, location = 0.0, alpha = 10.0, delta = 1.0, grid_n = 8001) {
  if (delta == 0.0) {
    return(rep(location, length(p)))
  }

  min_p <- max(min(p[p > 0.0], na.rm = TRUE), 1e-12)
  tail_width <- -log(min_p * 1e-4) / alpha
  L <- max(10.0 * sqrt(symmetric_hyperbolic_variance_from_z(alpha * delta, alpha)), delta + tail_width)
  y <- seq(-L, L, length.out = grid_n)
  density <- exp(-alpha * sqrt(delta^2 + y^2)) / (2.0 * delta * besselK(alpha * delta, 1))
  area <- c(0.0, cumsum(0.5 * (density[-1] + density[-length(density)]) * diff(y)))
  cdf <- area / area[length(area)]
  location + approx(cdf, y, xout = p, ties = "ordered", rule = 2)$y
}

hyperbolic_variance_from_delta <- function(delta, alpha, beta) {
  gamma <- sqrt(alpha^2 - beta^2)
  z <- delta * gamma
  k1 <- besselK(z, 1)
  k2 <- besselK(z, 2)
  k3 <- besselK(z, 3)
  delta / gamma * k2 / k1 +
    delta^2 * beta^2 / gamma^2 * (k3 / k1 - (k2 / k1)^2)
}

hyperbolic_solve_delta <- function(alpha, beta, v) {
  lo <- 1e-8
  min_v <- hyperbolic_variance_from_delta(lo, alpha, beta)
  if (v <= min_v) {
    stop("hyperbolic terminal_model cannot match sigma^2 * T with this hyperbolic_alpha and hyperbolic_beta; increase hyperbolic_alpha")
  }

  f <- function(delta) hyperbolic_variance_from_delta(delta, alpha, beta) - v
  hi <- max(1.0, sqrt(v))
  while (f(hi) < 0.0) {
    hi <- 2.0 * hi
  }
  uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

hyperbolic_log_mgf_centered <- function(t, alpha, beta, delta) {
  if (abs(beta + t) >= alpha) {
    return(Inf)
  }

  gamma0 <- sqrt(alpha^2 - beta^2)
  gamma_t <- sqrt(alpha^2 - (beta + t)^2)
  log(gamma0) - log(gamma_t) +
    log(besselK(delta * gamma_t, 1)) -
    log(besselK(delta * gamma0, 1))
}

hyperbolic_logreturn_params <- function(S0, mu, sigma, T, alpha = 10.0, beta = -2.0) {
  if (!is.finite(alpha) || !is.finite(beta) || alpha <= abs(beta)) {
    stop("hyperbolic terminal_model requires hyperbolic_alpha > abs(hyperbolic_beta)")
  }
  if (alpha <= abs(beta + 1.0)) {
    stop("hyperbolic terminal_model requires hyperbolic_alpha > abs(hyperbolic_beta + 1) for finite call prices")
  }

  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, alpha = alpha, beta = beta, delta = 0.0, v = v))
  }

  delta <- hyperbolic_solve_delta(alpha, beta, v)
  log_mgf_centered_at_one <- hyperbolic_log_mgf_centered(1.0, alpha, beta, delta)
  list(
    location = log(S0) + mu * T - log_mgf_centered_at_one,
    alpha = alpha,
    beta = beta,
    delta = delta,
    v = v
  )
}

qhyperbolic <- function(p, location = 0.0, alpha = 10.0, beta = -2.0, delta = 1.0, grid_n = 8001) {
  if (delta == 0.0) {
    return(rep(location, length(p)))
  }

  min_p <- max(min(p[p > 0.0], na.rm = TRUE), 1e-12)
  left_width <- -log(min_p * 1e-4) / (alpha + beta)
  right_width <- -log(min_p * 1e-4) / (alpha - beta)
  y <- seq(-max(delta + left_width, 1.0), max(delta + right_width, 1.0), length.out = grid_n)
  gamma <- sqrt(alpha^2 - beta^2)
  density <- gamma / (2.0 * alpha * delta * besselK(delta * gamma, 1)) *
    exp(-alpha * sqrt(delta^2 + y^2) + beta * y)
  area <- c(0.0, cumsum(0.5 * (density[-1] + density[-length(density)]) * diff(y)))
  cdf <- area / area[length(area)]
  location + approx(cdf, y, xout = p, ties = "ordered", rule = 2)$y
}

log_besselK <- function(x, nu) {
  log(besselK(x, nu, expon.scaled = TRUE)) - x
}

gh_log_mgf_centered <- function(t, lambda, alpha, beta, delta) {
  if (abs(beta + t) >= alpha) {
    return(Inf)
  }

  gamma0 <- sqrt(alpha^2 - beta^2)
  gamma_t <- sqrt(alpha^2 - (beta + t)^2)
  lambda * (log(gamma0) - log(gamma_t)) +
    log_besselK(delta * gamma_t, lambda) -
    log_besselK(delta * gamma0, lambda)
}

gh_variance_from_delta <- function(delta, lambda, alpha, beta) {
  domain <- alpha - abs(beta)
  h <- min(1e-4, 0.1 * domain)
  if (h <= 0.0) {
    return(Inf)
  }
  (gh_log_mgf_centered(h, lambda, alpha, beta, delta) -
     2.0 * gh_log_mgf_centered(0.0, lambda, alpha, beta, delta) +
     gh_log_mgf_centered(-h, lambda, alpha, beta, delta)) / h^2
}

gh_solve_delta <- function(lambda, alpha, beta, v, label = "generalized_hyperbolic") {
  lo <- 1e-6
  f <- function(delta) gh_variance_from_delta(delta, lambda, alpha, beta) - v
  if (f(lo) > 0.0) {
    stop(label, " terminal_model cannot match sigma^2 * T with these parameters; increase alpha or adjust beta/lambda")
  }
  hi <- max(1.0, sqrt(v))
  while (f(hi) < 0.0) {
    hi <- 2.0 * hi
  }
  uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

gh_logreturn_params <- function(S0, mu, sigma, T, lambda = 1.0, alpha = 10.0, beta = -2.0,
                                label = "generalized_hyperbolic") {
  if (!is.finite(lambda) || !is.finite(alpha) || !is.finite(beta) || alpha <= abs(beta)) {
    stop(label, " terminal_model requires finite lambda and alpha > abs(beta)")
  }
  if (alpha <= abs(beta + 1.0)) {
    stop(label, " terminal_model requires alpha > abs(beta + 1) for finite call prices")
  }

  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, lambda = lambda, alpha = alpha, beta = beta, delta = 0.0, v = v))
  }

  delta <- gh_solve_delta(lambda, alpha, beta, v, label = label)
  log_mgf_centered_at_one <- gh_log_mgf_centered(1.0, lambda, alpha, beta, delta)
  list(
    location = log(S0) + mu * T - log_mgf_centered_at_one,
    lambda = lambda,
    alpha = alpha,
    beta = beta,
    delta = delta,
    v = v
  )
}

gh_log_density_centered <- function(y, lambda, alpha, beta, delta) {
  gamma <- sqrt(alpha^2 - beta^2)
  r <- sqrt(delta^2 + y^2)
  lambda * (log(gamma) - log(delta)) -
    0.5 * log(2.0 * pi) -
    log_besselK(delta * gamma, lambda) +
    log_besselK(alpha * r, lambda - 0.5) +
    beta * y -
    (0.5 - lambda) * (log(r) - log(alpha))
}

qgh <- function(p, location = 0.0, lambda = 1.0, alpha = 10.0, beta = -2.0, delta = 1.0, grid_n = 8001) {
  if (delta == 0.0) {
    return(rep(location, length(p)))
  }

  min_p <- max(min(p[p > 0.0], na.rm = TRUE), 1e-12)
  left_width <- -log(min_p * 1e-4) / (alpha + beta)
  right_width <- -log(min_p * 1e-4) / (alpha - beta)
  y <- seq(-max(delta + left_width, 1.0), max(delta + right_width, 1.0), length.out = grid_n)
  log_density <- gh_log_density_centered(y, lambda, alpha, beta, delta)
  density <- exp(log_density - max(log_density, na.rm = TRUE))
  area <- c(0.0, cumsum(0.5 * (density[-1] + density[-length(density)]) * diff(y)))
  cdf <- area / area[length(area)]
  location + approx(cdf, y, xout = p, ties = "ordered", rule = 2)$y
}

vg_variance_from_alpha <- function(alpha, beta, shape) {
  2.0 * shape * (alpha^2 + beta^2) / (alpha^2 - beta^2)^2
}

vg_solve_alpha <- function(beta, shape, v) {
  lo <- abs(beta) + 1.0 + 1e-8
  f <- function(alpha) vg_variance_from_alpha(alpha, beta, shape) - v
  hi <- max(lo * 2.0, 2.0)
  while (f(hi) > 0.0) {
    hi <- 2.0 * hi
  }
  uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

vg_logreturn_params <- function(S0, mu, sigma, T, shape = 1.0, beta = -0.2) {
  if (!is.finite(shape) || shape <= 0.0 || !is.finite(beta)) {
    stop("variance_gamma terminal_model requires vg_shape > 0 and finite vg_beta")
  }

  v <- sigma^2 * T
  if (v <= 0.0) {
    return(list(location = log(S0) + mu * T, shape = shape, alpha = abs(beta) + 2.0, beta = beta, v = v))
  }

  alpha <- vg_solve_alpha(beta, shape, v)
  log_mgf_centered_at_one <- shape * (log(alpha^2 - beta^2) - log(alpha^2 - (beta + 1.0)^2))
  list(
    location = log(S0) + mu * T - log_mgf_centered_at_one,
    shape = shape,
    alpha = alpha,
    beta = beta,
    v = v
  )
}

vg_log_density_centered <- function(y, shape, alpha, beta) {
  ay <- pmax(abs(y), 1e-12)
  shape * log(alpha^2 - beta^2) +
    (shape - 0.5) * log(ay) +
    log_besselK(alpha * ay, shape - 0.5) +
    beta * y -
    0.5 * log(pi) -
    lgamma(shape) -
    (shape - 0.5) * log(2.0 * alpha)
}

qvg <- function(p, location = 0.0, shape = 1.0, alpha = 10.0, beta = -0.2, grid_n = 8001) {
  min_p <- max(min(p[p > 0.0], na.rm = TRUE), 1e-12)
  left_width <- -log(min_p * 1e-4) / (alpha + beta)
  right_width <- -log(min_p * 1e-4) / (alpha - beta)
  y <- seq(-max(left_width, 1.0), max(right_width, 1.0), length.out = grid_n)
  log_density <- vg_log_density_centered(y, shape, alpha, beta)
  density <- exp(log_density - max(log_density, na.rm = TRUE))
  area <- c(0.0, cumsum(0.5 * (density[-1] + density[-length(density)]) * diff(y)))
  cdf <- area / area[length(area)]
  location + approx(cdf, y, xout = p, ties = "ordered", rule = 2)$y
}

lognormal_mixture_cdf <- function(x, weights, m, v) {
  if (x <= 0.0) {
    return(0.0)
  }
  sum(weights * pnorm((log(x) - m) / sqrt(v)))
}

lognormal_mixture_quantile <- function(p, weights, m, v) {
  if (p <= 0.0) {
    return(0.0)
  }
  if (p >= 1.0) {
    return(Inf)
  }

  component_q <- exp(m + sqrt(v) * qnorm(pmax(p, 1e-12)))
  hi <- max(component_q, 1.0)
  while (lognormal_mixture_cdf(hi, weights, m, v) < p) {
    hi <- 2.0 * hi
  }

  uniroot(
    function(x) lognormal_mixture_cdf(x, weights, m, v) - p,
    lower = 0.0,
    upper = hi,
    tol = 1e-10
  )$root
}

terminal_price_quantiles <- function(model,
                                     S0,
                                     mu,
                                     sigma,
                                     T,
                                     n,
                                     terminal_floor = NA_real_,
                                     mixture_weights = NULL,
                                     mixture_mu = NULL,
                                     mixture_sigma = NULL,
                                     hyperbolic_alpha = 10.0,
                                     hyperbolic_beta = -2.0,
                                     gh_lambda = 1.0,
                                     gh_alpha = 10.0,
                                     gh_beta = -2.0,
                                     nig_alpha = 10.0,
                                     nig_beta = -2.0,
                                     vg_shape = 1.0,
                                     vg_beta = -0.2,
                                     probs = NULL) {
  if (is.null(probs)) {
    probs <- seq(0.001, 0.999, length.out = n)
  }

  if (model == "lognormal") {
    par <- ln_params_from_gbm(S0, mu, sigma, T)
    S <- exp(par$m + sqrt(par$v) * qnorm(probs))
  } else if (model == "logistic") {
    par <- logistic_logreturn_params(S0, mu, sigma, T)
    if (par$scale == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qlogis(probs, location = par$location, scale = par$scale))
    }
  } else if (model == "hyperbolic_secant") {
    par <- hyperbolic_secant_logreturn_params(S0, mu, sigma, T)
    if (par$scale == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qhyperbolic_secant(probs, location = par$location, scale = par$scale))
    }
  } else if (model == "symmetric_hyperbolic") {
    par <- symmetric_hyperbolic_logreturn_params(S0, mu, sigma, T, alpha = hyperbolic_alpha)
    if (par$delta == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qsymmetric_hyperbolic(probs, location = par$location, alpha = par$alpha, delta = par$delta))
    }
  } else if (model == "hyperbolic") {
    par <- hyperbolic_logreturn_params(S0, mu, sigma, T, alpha = hyperbolic_alpha, beta = hyperbolic_beta)
    if (par$delta == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qhyperbolic(probs, location = par$location, alpha = par$alpha, beta = par$beta, delta = par$delta))
    }
  } else if (model == "generalized_hyperbolic") {
    par <- gh_logreturn_params(S0, mu, sigma, T, lambda = gh_lambda, alpha = gh_alpha, beta = gh_beta)
    if (par$delta == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qgh(probs, location = par$location, lambda = par$lambda, alpha = par$alpha, beta = par$beta, delta = par$delta))
    }
  } else if (model == "normal_inverse_gaussian") {
    par <- gh_logreturn_params(S0, mu, sigma, T, lambda = -0.5, alpha = nig_alpha, beta = nig_beta,
                               label = "normal_inverse_gaussian")
    if (par$delta == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qgh(probs, location = par$location, lambda = par$lambda, alpha = par$alpha, beta = par$beta, delta = par$delta))
    }
  } else if (model == "variance_gamma") {
    par <- vg_logreturn_params(S0, mu, sigma, T, shape = vg_shape, beta = vg_beta)
    if (par$v == 0.0) {
      S <- rep(exp(par$location), length(probs))
    } else {
      S <- exp(qvg(probs, location = par$location, shape = par$shape, alpha = par$alpha, beta = par$beta))
    }
  } else if (model == "normal") {
    terminal_mean <- S0 * exp(mu * T)
    terminal_sd <- S0 * sigma * sqrt(T)
    S <- qnorm(probs, mean = terminal_mean, sd = terminal_sd)
  } else if (model == "lognormal_mixture") {
    if (is.null(mixture_weights) || is.null(mixture_mu) || is.null(mixture_sigma)) {
      stop("lognormal_mixture requires mixture_weights, mixture_mu, and mixture_sigma")
    }
    stopifnot(length(mixture_weights) == length(mixture_mu))
    stopifnot(length(mixture_weights) == length(mixture_sigma))
    weights <- normalize_mixture_weights(mixture_weights)
    par <- lognormal_mixture_params(S0, mixture_mu, mixture_sigma, T)
    S <- vapply(
      probs,
      function(p) lognormal_mixture_quantile(p, weights, par$m, par$v),
      numeric(1)
    )
  } else {
    stop("terminal_model must be 'lognormal', 'logistic', 'hyperbolic_secant', 'symmetric_hyperbolic', 'hyperbolic', 'generalized_hyperbolic', 'normal_inverse_gaussian', 'variance_gamma', 'normal', or 'lognormal_mixture'")
  }

  apply_terminal_floor(S, terminal_floor)
}

# E[((S_T - K1)^+) * ((S_T - K2)^+)]
call_payoff_cross_moment <- function(m, v, K1, K2) {
  Kstar <- max(K1, K2)

  m2 <- ln_trunc_moment(m, v, Kstar, 2)
  m1 <- ln_trunc_moment(m, v, Kstar, 1)
  p0 <- Phi(ln_d(m, v, Kstar, 0))

  m2 - (K1 + K2) * m1 + K1 * K2 * p0
}

# Cov[((S_T - K1)^+), ((S_T - K2)^+)]
call_payoff_covariance <- function(m, v, K1, K2) {
  e12 <- call_payoff_cross_moment(m, v, K1, K2)
  e1 <- call_payoff_mean(m, v, K1)
  e2 <- call_payoff_mean(m, v, K2)

  e12 - e1 * e2
}

# Mean vector for calls with strikes K
call_payoff_mean_vec <- function(m, v, K) {
  sapply(K, function(k) call_payoff_mean(m, v, k))
}

# Covariance matrix for calls with strikes K
call_payoff_cov_mat <- function(m, v, K) {
  n <- length(K)
  out <- matrix(0.0, n, n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      out[i, j] <- call_payoff_covariance(m, v, K[i], K[j])
    }
  }

  rownames(out) <- paste0("K=", K)
  colnames(out) <- paste0("K=", K)
  out
}

# Correlation matrix from covariance matrix
call_payoff_cor_mat <- function(covmat) {
  sd <- sqrt(diag(covmat))
  covmat / outer(sd, sd)
}

# Black-Scholes call price under continuous dividend yield q
bs_call_price <- function(S0, K, r, q, sigma, T) {
  if (K <= 0) {
    return(S0 * exp(-q * T))
  }

  vol_sqrt_t <- sigma * sqrt(T)
  d1 <- (log(S0 / K) + (r - q + 0.5 * sigma^2) * T) / vol_sqrt_t
  d2 <- d1 - vol_sqrt_t
  S0 * exp(-q * T) * Phi(d1) - K * exp(-r * T) * Phi(d2)
}

bs_call_price_vec <- function(S0, K, r, q, sigma, T) {
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(K))
  }

  stopifnot(length(K) == length(sigma))
  sapply(seq_along(K), function(i) bs_call_price(S0, K[i], r, q, sigma[i], T))
}

bs_put_price <- function(S0, K, r, q, sigma, T) {
  bs_call_price(S0, K, r, q, sigma, T) -
    S0 * exp(-q * T) +
    K * exp(-r * T)
}

bs_option_price_vec <- function(S0, K, r, q, sigma, T, type) {
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(K))
  }

  stopifnot(length(K) == length(sigma))
  stopifnot(length(K) == length(type))

  sapply(seq_along(K), function(i) {
    if (type[i] == "call") {
      bs_call_price(S0, K[i], r, q, sigma[i], T)
    } else if (type[i] == "put") {
      bs_put_price(S0, K[i], r, q, sigma[i], T)
    } else {
      stop("type must be 'call' or 'put'")
    }
  })
}

bs_option_greeks <- function(S0, K, r, q, sigma, T, type) {
  if (K <= 0 && type == "call") {
    return(c(delta = exp(-q * T), gamma = 0.0, vega = 0.0))
  }
  if (!is.finite(sigma) || sigma <= 0.0 || T <= 0.0 || K <= 0.0) {
    return(c(delta = NA_real_, gamma = NA_real_, vega = NA_real_))
  }

  vol_sqrt_t <- sigma * sqrt(T)
  d1 <- (log(S0 / K) + (r - q + 0.5 * sigma^2) * T) / vol_sqrt_t
  discounted_spot <- exp(-q * T)
  delta <- if (type == "call") {
    discounted_spot * Phi(d1)
  } else if (type == "put") {
    discounted_spot * (Phi(d1) - 1.0)
  } else {
    stop("type must be 'call' or 'put'")
  }
  gamma <- discounted_spot * dnorm(d1) / (S0 * vol_sqrt_t)
  vega <- S0 * discounted_spot * dnorm(d1) * sqrt(T)
  c(delta = delta, gamma = gamma, vega = vega)
}

bs_option_greeks_table <- function(S0, K, r, q, sigma, T, type) {
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(K))
  }
  stopifnot(length(K) == length(sigma))
  stopifnot(length(K) == length(type))

  out <- t(vapply(seq_along(K), function(i) {
    bs_option_greeks(S0, K[i], r, q, sigma[i], T, type[i])
  }, numeric(3)))
  out <- as.data.frame(out)
  rownames(out) <- NULL
  out
}

var_spot_scenarios <- function(S0,
                               mu,
                               sigma,
                               horizon_days = 1.0,
                               trading_days = 252.0,
                               n_scenarios = 1001,
                               return_model = "normal_log",
                               terminal_floor = 0.0) {
  dt <- horizon_days / trading_days
  probs <- (seq_len(n_scenarios) - 0.5) / n_scenarios
  z <- qnorm(probs)

  if (return_model == "normal_log") {
    log_return <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * z
    S <- S0 * exp(log_return)
  } else if (return_model == "normal_simple") {
    simple_return <- mu * dt + sigma * sqrt(dt) * z
    S <- S0 * (1.0 + simple_return)
  } else {
    stop("var_return_model must be 'normal_log' or 'normal_simple'")
  }

  apply_terminal_floor(S, terminal_floor)
}

option_mtm_price_matrix <- function(S,
                                    K,
                                    r,
                                    q,
                                    sigma,
                                    T,
                                    type) {
  if (T <= 0.0) {
    return(option_payoff_matrix(S, K, type))
  }

  out <- matrix(0.0, nrow = length(S), ncol = length(K))
  for (i in seq_along(S)) {
    out[i, ] <- bs_option_price_vec(S[i], K, r, q, sigma, T, type)
  }
  colnames(out) <- paste(type, K, sep = "_")
  out
}

portfolio_mtm_var_stats <- function(contracts,
                                    cash,
                                    current_mid_price,
                                    S0,
                                    K,
                                    r,
                                    q,
                                    mid_vol,
                                    T,
                                    type,
                                    var_mu,
                                    var_sigma,
                                    horizon_days = 1.0,
                                    trading_days = 252.0,
                                    n_scenarios = 1001,
                                    return_model = "normal_log",
                                    terminal_floor = 0.0,
                                    conf_levels = c(0.95, 0.99)) {
  dt <- horizon_days / trading_days
  S <- var_spot_scenarios(
    S0 = S0,
    mu = var_mu,
    sigma = var_sigma,
    horizon_days = horizon_days,
    trading_days = trading_days,
    n_scenarios = n_scenarios,
    return_model = return_model,
    terminal_floor = terminal_floor
  )
  future_T <- max(T - dt, 0.0)
  future_prices <- option_mtm_price_matrix(
    S = S,
    K = K,
    r = r,
    q = q,
    sigma = mid_vol,
    T = future_T,
    type = type
  )
  current_value <- cash + sum(current_mid_price * contracts)
  future_value <- cash * exp(r * dt) + as.numeric(future_prices %*% contracts)
  pnl <- future_value - current_value
  loss <- -pnl

  out <- list(
    current_mtm_wealth = current_value,
    mean_mtm_pnl = mean(pnl),
    sd_mtm_pnl = sd(pnl)
  )
  for (level in conf_levels) {
    suffix <- if (abs(100.0 * level - round(100.0 * level)) < 1e-8) {
      as.character(as.integer(round(100.0 * level)))
    } else {
      gsub("\\.", "p", format(100.0 * level, trim = TRUE, scientific = FALSE))
    }
    threshold <- as.numeric(quantile(loss, level, names = FALSE, type = 8))
    tail_loss <- loss[loss >= threshold]
    out[[paste0("var_", suffix)]] <- threshold
    out[[paste0("es_", suffix)]] <- mean(tail_loss)
  }
  out
}

multistock_var_price_scenarios <- function(S0,
                                           mu,
                                           cov_matrix,
                                           horizon_days = 1.0,
                                           trading_days = 252.0,
                                           n_scenarios = 1001,
                                           return_model = "normal_log",
                                           terminal_floor = 0.0,
                                           seed = NA_integer_) {
  S0 <- as.numeric(S0)
  mu <- as.numeric(mu)
  stopifnot(length(S0) == length(mu))
  stopifnot(nrow(cov_matrix) == length(S0))
  stopifnot(ncol(cov_matrix) == length(S0))
  if (!is.na(seed)) {
    set.seed(as.integer(seed))
  }
  dt <- horizon_days / trading_days
  n_assets <- length(S0)
  chol_cov <- chol(cov_matrix * dt)
  z <- matrix(rnorm(n_scenarios * n_assets), nrow = n_scenarios, ncol = n_assets)
  shocks <- z %*% chol_cov

  if (return_model == "normal_log") {
    variance <- diag(cov_matrix)
    drift <- matrix(
      log(S0) + (mu - 0.5 * variance) * dt,
      nrow = n_scenarios,
      ncol = n_assets,
      byrow = TRUE
    )
    S <- exp(drift + shocks)
  } else if (return_model == "normal_simple") {
    drift <- matrix(mu * dt, nrow = n_scenarios, ncol = n_assets, byrow = TRUE)
    S <- matrix(S0, nrow = n_scenarios, ncol = n_assets, byrow = TRUE) * (1.0 + drift + shocks)
  } else {
    stop("var_return_model must be 'normal_log' or 'normal_simple'")
  }

  apply_terminal_floor(S, terminal_floor)
}

multistock_mtm_price_matrix <- function(spot_scenarios,
                                        underlying,
                                        K,
                                        r,
                                        q,
                                        sigma,
                                        T,
                                        type,
                                        underlying_names) {
  if (T <= 0.0) {
    return(multistock_payoff_matrix(spot_scenarios, underlying, K, type, underlying_names))
  }

  out <- matrix(0.0, nrow = nrow(spot_scenarios), ncol = length(K))
  for (j in seq_along(K)) {
    asset_idx <- match(underlying[j], underlying_names)
    if (is.na(asset_idx)) {
      stop("Unknown underlying: ", underlying[j])
    }
    S <- spot_scenarios[, asset_idx]
    if (K[j] <= 0.0 && type[j] == "call") {
      out[, j] <- S
    } else if (type[j] == "call") {
      out[, j] <- bs_call_price(S, K[j], r, q[asset_idx], sigma[j], T)
    } else if (type[j] == "put") {
      out[, j] <- bs_put_price(S, K[j], r, q[asset_idx], sigma[j], T)
    } else {
      stop("type must be 'call' or 'put'")
    }
  }
  colnames(out) <- paste(underlying, type, K, sep = "_")
  out
}

multistock_portfolio_mtm_var_stats <- function(contracts,
                                               cash,
                                               current_mid_price,
                                               S0,
                                               underlying,
                                               K,
                                               r,
                                               q,
                                               mid_vol,
                                               T,
                                               type,
                                               var_mu,
                                               var_cov_matrix,
                                               underlying_names,
                                               horizon_days = 1.0,
                                               trading_days = 252.0,
                                               n_scenarios = 1001,
                                               return_model = "normal_log",
                                               terminal_floor = 0.0,
                                               conf_levels = c(0.95, 0.99),
                                               seed = NA_integer_) {
  dt <- horizon_days / trading_days
  spot_scenarios <- multistock_var_price_scenarios(
    S0 = S0,
    mu = var_mu,
    cov_matrix = var_cov_matrix,
    horizon_days = horizon_days,
    trading_days = trading_days,
    n_scenarios = n_scenarios,
    return_model = return_model,
    terminal_floor = terminal_floor,
    seed = seed
  )
  colnames(spot_scenarios) <- underlying_names
  future_T <- max(T - dt, 0.0)
  future_prices <- multistock_mtm_price_matrix(
    spot_scenarios = spot_scenarios,
    underlying = underlying,
    K = K,
    r = r,
    q = q,
    sigma = mid_vol,
    T = future_T,
    type = type,
    underlying_names = underlying_names
  )
  current_value <- cash + sum(current_mid_price * contracts)
  future_value <- cash * exp(r * dt) + as.numeric(future_prices %*% contracts)
  pnl <- future_value - current_value
  loss <- -pnl

  out <- list(
    current_mtm_wealth = current_value,
    mean_mtm_pnl = mean(pnl),
    sd_mtm_pnl = sd(pnl)
  )
  for (level in conf_levels) {
    suffix <- if (abs(100.0 * level - round(100.0 * level)) < 1e-8) {
      as.character(as.integer(round(100.0 * level)))
    } else {
      gsub("\\.", "p", format(100.0 * level, trim = TRUE, scientific = FALSE))
    }
    threshold <- as.numeric(quantile(loss, level, names = FALSE, type = 8))
    tail_loss <- loss[loss >= threshold]
    out[[paste0("var_", suffix)]] <- threshold
    out[[paste0("es_", suffix)]] <- mean(tail_loss)
  }
  out
}

mtm_constraint_col_names <- function(prefix, conf_levels) {
  pct <- 100 * conf_levels
  pct_label <- ifelse(abs(pct - round(pct)) < 1e-8,
                      as.character(as.integer(round(pct))),
                      gsub("\\.", "p", format(pct, trim = TRUE, scientific = FALSE)))
  paste0(prefix, "_", pct_label)
}

normalize_mtm_loss_constraints <- function(x, label) {
  if (is.null(x)) {
    return(data.frame(conf_level = numeric(), max_loss = numeric()))
  }
  if (!is.data.frame(x)) {
    stop(label, " must be a data.frame with conf_level and max_loss")
  }
  if (nrow(x) == 0) {
    return(data.frame(conf_level = numeric(), max_loss = numeric()))
  }
  if (!all(c("conf_level", "max_loss") %in% names(x))) {
    stop(label, " must contain conf_level and max_loss")
  }
  out <- data.frame(
    conf_level = as.numeric(x$conf_level),
    max_loss = as.numeric(x$max_loss)
  )
  if (any(!is.finite(out$conf_level) | out$conf_level <= 0.0 | out$conf_level >= 1.0)) {
    stop(label, " conf_level values must be in (0, 1)")
  }
  if (any(!is.finite(out$max_loss))) {
    stop(label, " max_loss values must be finite")
  }
  out
}

simulate_gbm_paths <- function(S0,
                               mu,
                               sigma,
                               T,
                               steps,
                               n_paths,
                               seed = NA_integer_) {
  if (steps < 1L) {
    stop("steps must be positive")
  }
  if (!is.na(seed)) {
    set.seed(as.integer(seed))
  }
  dt <- T / steps
  z <- matrix(rnorm(n_paths * steps), nrow = n_paths, ncol = steps)
  log_returns <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * z
  log_paths <- cbind(0.0, t(apply(log_returns, 1L, cumsum)))
  S0 * exp(log_paths)
}

portfolio_delta_at_state <- function(S,
                                     K,
                                     r,
                                     q,
                                     sigma,
                                     remaining_T,
                                     type,
                                     contracts) {
  if (remaining_T <= 0.0) {
    return(0.0)
  }
  greeks <- bs_option_greeks_table(S, K, r, q, sigma, remaining_T, type)
  sum(contracts * greeks$delta, na.rm = TRUE)
}

delta_hedged_portfolio_stats <- function(contracts,
                                         cash,
                                         budget,
                                         S0,
                                         K,
                                         r,
                                         q,
                                         mid_vol,
                                         T,
                                         type,
                                         hedge_steps = 0L,
                                         hedge_paths = 1000L,
                                         hedge_mu = r - q,
                                         hedge_sigma,
                                         hedge_seed = NA_integer_,
                                         stock_transaction_cost = 0.0) {
  hedge_steps <- as.integer(hedge_steps)
  hedge_paths <- as.integer(hedge_paths)
  if (hedge_steps <= 0L) {
    return(list(
      hedged_mean_wealth = NA_real_,
      hedged_sd_wealth = NA_real_,
      hedged_sharpe = NA_real_,
      hedged_skew = NA_real_,
      hedged_ex.kurt = NA_real_,
      hedged_mean_pnl = NA_real_,
      hedged_sd_pnl = NA_real_,
      hedged_avg_abs_final_delta = NA_real_
    ))
  }
  if (hedge_paths <= 1L) {
    stop("delta_hedge_paths must be greater than 1")
  }
  if (!is.finite(hedge_sigma) || hedge_sigma < 0.0) {
    stop("delta_hedge_sigma must be finite and nonnegative")
  }

  paths <- simulate_gbm_paths(
    S0 = S0,
    mu = hedge_mu,
    sigma = hedge_sigma,
    T = T,
    steps = hedge_steps + 1L,
    n_paths = hedge_paths,
    seed = hedge_seed
  )
  dt <- T / (hedge_steps + 1L)
  cash_path <- rep(cash, hedge_paths)
  hedge_shares <- numeric(hedge_paths)

  for (step in 0:hedge_steps) {
    S_step <- paths[, step + 1L]
    remaining_T <- T - step * dt
    target_hedge <- -vapply(S_step, function(S) {
      portfolio_delta_at_state(
        S = S,
        K = K,
        r = r,
        q = q,
        sigma = mid_vol,
        remaining_T = remaining_T,
        type = type,
        contracts = contracts
      )
    }, numeric(1))
    trade_shares <- target_hedge - hedge_shares
    cash_path <- cash_path - trade_shares * S_step -
      abs(trade_shares) * stock_transaction_cost
    hedge_shares <- target_hedge
    cash_path <- cash_path * exp(r * dt)
  }

  S_T <- paths[, ncol(paths)]
  payoff <- as.numeric(option_payoff_matrix(S_T, K, type) %*% contracts)
  terminal_wealth <- cash_path + hedge_shares * S_T + payoff
  rf_terminal <- budget * exp(r * T)
  pnl <- terminal_wealth - rf_terminal
  sd_wealth <- sd(terminal_wealth)
  mean_wealth <- mean(terminal_wealth)
  centered <- terminal_wealth - mean_wealth
  skew <- if (sd_wealth > 0.0) mean(centered^3) / sd_wealth^3 else NaN
  ex_kurt <- if (sd_wealth > 0.0) mean(centered^4) / sd_wealth^4 - 3.0 else NaN
  final_delta <- vapply(S_T, function(S) {
    portfolio_delta_at_state(
      S = S,
      K = K,
      r = r,
      q = q,
      sigma = mid_vol,
      remaining_T = 0.0,
      type = type,
      contracts = contracts
    )
  }, numeric(1)) + hedge_shares

  list(
    hedged_mean_wealth = mean_wealth,
    hedged_sd_wealth = sd_wealth,
    hedged_sharpe = (mean_wealth - rf_terminal) / sd_wealth,
    hedged_skew = skew,
    hedged_ex.kurt = ex_kurt,
    hedged_mean_pnl = mean(pnl),
    hedged_sd_pnl = sd(pnl),
    hedged_avg_abs_final_delta = mean(abs(final_delta))
  )
}

option_holding_period_sharpe <- function(S0,
                                         K,
                                         r,
                                         q,
                                         realized_sigma,
                                         T,
                                         type = "call",
                                         price = NA_real_,
                                         implied_sigma = NA_real_,
                                         mu = r - q) {
  if (any(type == "straddle")) {
    if (!all(type == "straddle")) {
      stop("For ann_sharpe(), do not mix type='straddle' with call/put rows in one call")
    }
    if (is.na(price[1])) {
      if (is.na(implied_sigma[1])) {
        stop("Provide either price or implied_sigma")
      }
      call_price <- bs_option_price_vec(S0, K, r, q, implied_sigma, T, rep("call", length(K)))
      put_price <- bs_option_price_vec(S0, K, r, q, implied_sigma, T, rep("put", length(K)))
      price <- call_price + put_price
    }

    par <- ln_params_from_gbm(S0, mu, realized_sigma, T)
    expanded_K <- rep(K, each = 2)
    expanded_type <- rep(c("call", "put"), times = length(K))
    payoff_stats <- vapply(seq_along(K), function(i) {
      stats <- option_portfolio_stats(
        par$m,
        par$v,
        expanded_K[(2 * i - 1):(2 * i)],
        expanded_type[(2 * i - 1):(2 * i)],
        c(1.0, 1.0)
      )
      stats[c("mean", "sd", "skew", "ex.kurt")]
    }, numeric(4))
    payoff_stats <- t(payoff_stats)
    holding_sharpe <- (payoff_stats[, "mean"] - price * exp(r * T)) / payoff_stats[, "sd"]

    out <- data.frame(
      type = type,
      strike = K,
      price = as.numeric(price),
      implied_sigma = implied_sigma,
      mu = mu,
      realized_sigma = realized_sigma,
      T = T,
      mean_payoff = payoff_stats[, "mean"],
      sd_payoff = payoff_stats[, "sd"],
      mean_profit = payoff_stats[, "mean"] - price * exp(r * T),
      sd_profit = payoff_stats[, "sd"],
      sharpe = as.numeric(holding_sharpe)
    )
    rownames(out) <- NULL
    return(out)
  }

  if (is.na(price)) {
    if (is.na(implied_sigma)) {
      stop("Provide either price or implied_sigma")
    }
    price <- bs_option_price_vec(S0, K, r, q, implied_sigma, T, type)
  }

  par <- ln_params_from_gbm(S0, mu, realized_sigma, T)
  stats <- option_payoff_stats_table(par$m, par$v, K, type)
  holding_sharpe <- (stats$mean - price * exp(r * T)) / stats$sd

  data.frame(
    type = type,
    strike = K,
    price = as.numeric(price),
    implied_sigma = implied_sigma,
    mu = mu,
    realized_sigma = realized_sigma,
    T = T,
    mean_payoff = stats$mean,
    sd_payoff = stats$sd,
    mean_profit = stats$mean - price * exp(r * T),
    sd_profit = stats$sd,
    sharpe = as.numeric(holding_sharpe)
  )
}

option_annualized_sharpe <- function(S0,
                                     K,
                                     r,
                                     q = 0.0,
                                     realized_sigma,
                                     T,
                                     type = "call",
                                     price = NA_real_,
                                     implied_sigma = NA_real_,
                                     mu = r - q) {
  out <- option_holding_period_sharpe(
    S0 = S0,
    K = K,
    r = r,
    q = q,
    realized_sigma = realized_sigma,
    T = T,
    type = type,
    price = price,
    implied_sigma = implied_sigma,
    mu = mu
  )
  out$ann_sharpe <- out$sharpe / sqrt(T)
  out
}

ann_sharpe <- option_annualized_sharpe

bs_option_implied_vol <- function(target_price,
                                  S0,
                                  K,
                                  r,
                                  q,
                                  T,
                                  type,
                                  lower = 1e-6,
                                  upper = 5.0,
                                  tol = 1e-8) {
  if (!is.finite(target_price) || K <= 0) {
    return(NA_real_)
  }

  f <- function(sig) {
    bs_option_price_vec(S0, K, r, q, sig, T, type) - target_price
  }
  flo <- f(lower)
  fhi <- f(upper)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0.0) {
    return(NA_real_)
  }

  uniroot(f, lower = lower, upper = upper, tol = tol)$root
}

bs_option_implied_vol_vec <- function(target_price, S0, K, r, q, T, type) {
  stopifnot(length(target_price) == length(K))
  stopifnot(length(K) == length(type))
  vapply(
    seq_along(K),
    function(i) bs_option_implied_vol(target_price[i], S0, K[i], r, q, T, type[i]),
    numeric(1)
  )
}

vol_bid_ask_from_quadratic_spread <- function(S0,
                                              K,
                                              T,
                                              r,
                                              q,
                                              mid_vol,
                                              spread_atm,
                                              spread_quad,
                                              min_bid_vol = 1e-4) {
  forward <- S0 * exp((r - q) * T)
  x <- rep(0.0, length(K))
  option_idx <- K > 0
  x[option_idx] <- log(K[option_idx] / forward)

  vol_spread <- rep(0.0, length(K))
  vol_spread[option_idx] <- spread_atm + spread_quad * x[option_idx]^2
  bid_vol <- mid_vol
  ask_vol <- mid_vol
  bid_vol[option_idx] <- pmax(mid_vol[option_idx] - 0.5 * vol_spread[option_idx], min_bid_vol)
  ask_vol[option_idx] <- mid_vol[option_idx] + 0.5 * vol_spread[option_idx]

  data.frame(
    mid_vol = mid_vol,
    bid_vol = bid_vol,
    ask_vol = ask_vol,
    vol_spread = vol_spread
  )
}

position_exec_prices <- function(contracts, bid_prices, ask_prices) {
  ifelse(contracts < 0.0, bid_prices, ask_prices)
}

executable_trade_cost <- function(contracts, bid_prices, ask_prices) {
  sum(position_exec_prices(contracts, bid_prices, ask_prices) * contracts)
}

rebalance_trade_contracts <- function(final_contracts, initial_contracts) {
  final_contracts - initial_contracts
}

rebalance_trade_cost <- function(final_contracts, initial_contracts, bid_prices, ask_prices) {
  executable_trade_cost(
    rebalance_trade_contracts(final_contracts, initial_contracts),
    bid_prices,
    ask_prices
  )
}

check_rebalance_trade_costs <- function(tol = 1e-10) {
  bid <- 9.0
  ask <- 11.0
  checks <- c(
    no_trade = rebalance_trade_cost(1.0, 1.0, bid, ask) - 0.0,
    buy_more = rebalance_trade_cost(1.5, 1.0, bid, ask) - 0.5 * ask,
    sell_down = rebalance_trade_cost(0.25, 1.0, bid, ask) - (-0.75 * bid),
    cover_short = rebalance_trade_cost(-0.25, -1.0, bid, ask) - 0.75 * ask,
    short_more = rebalance_trade_cost(-1.5, -1.0, bid, ask) - (-0.5 * bid)
  )

  if (any(abs(checks) > tol)) {
    stop("rebalance trade cost sanity check failed")
  }

  TRUE
}

modified_sharpe_ratio <- function(sharpe, skew, ex.kurt) {
  sharpe * (1 + skew * sharpe / 6 - ex.kurt * sharpe^2 / 24)
}

option_payoff_dependency_matrix <- function(payoff_scenarios, tol = 1e-12) {
  positive_payoff <- payoff_scenarios > 0.0
  payout_prob <- colMeans(positive_payoff)
  joint_prob <- crossprod(positive_payoff) / nrow(positive_payoff)
  n <- length(payout_prob)
  lambda <- matrix(0.0, n, n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      denom <- payout_prob[i] * payout_prob[j]
      lambda[i, j] <- if (denom > tol) joint_prob[i, j] / denom else 0.0
    }
  }

  lambda <- 0.5 * (lambda + t(lambda))
  rownames(lambda) <- colnames(payoff_scenarios)
  colnames(lambda) <- colnames(payoff_scenarios)
  lambda
}

psd_repair_matrix <- function(x, eigen_floor = 1e-8) {
  x <- 0.5 * (x + t(x))
  eig <- eigen(x, symmetric = TRUE)
  raw_values <- eig$values
  repaired_values <- pmax(raw_values, eigen_floor)
  repaired <- eig$vectors %*% diag(repaired_values, nrow = length(repaired_values)) %*%
    t(eig$vectors)
  repaired <- 0.5 * (repaired + t(repaired))
  rownames(repaired) <- rownames(x)
  colnames(repaired) <- colnames(x)

  list(
    matrix = repaired,
    raw_min_eigen = min(raw_values),
    raw_max_eigen = max(raw_values),
    repaired_min_eigen = min(repaired_values),
    repaired_max_eigen = max(repaired_values),
    condition_number = max(repaired_values) / max(min(repaired_values), .Machine$double.eps),
    repaired = any(raw_values < eigen_floor)
  )
}

long_only_mean_variance_opt <- function(edge,
                                        covmat,
                                        risk_aversion,
                                        max_weight = 1.0,
                                        force_full_investment = FALSE,
                                        tol = 1e-10) {
  n <- length(edge)
  best_w <- numeric(n)
  best_obj <- if (force_full_investment) -Inf else 0.0

  objective <- function(w) {
    sum(edge * w) - risk_aversion * as.numeric(t(w) %*% covmat %*% w)
  }

  try_candidate <- function(idx, w_active) {
    if (any(!is.finite(w_active))) {
      return()
    }
    if (any(w_active < -tol)) {
      return()
    }

    w <- numeric(n)
    w[idx] <- pmax(w_active, 0.0)

    if (force_full_investment) {
      s <- sum(w)
      if (s <= tol) {
        return()
      }
      w <- max_weight * w / s
    } else if (sum(w) > max_weight + 1e-8) {
      return()
    }

    obj <- objective(w)
    if (is.finite(obj) && obj > best_obj + 1e-12) {
      best_w <<- w
      best_obj <<- obj
    }
  }

  for (active_size in seq_len(n)) {
    active_sets <- combn(seq_len(n), active_size, simplify = FALSE)

    for (idx in active_sets) {
      sigma <- covmat[idx, idx, drop = FALSE]
      e <- edge[idx]

      if (!force_full_investment) {
        w_uncon <- try(qr.solve(2.0 * risk_aversion * sigma, e), silent = TRUE)
        if (!inherits(w_uncon, "try-error")) {
          try_candidate(idx, as.numeric(w_uncon))
        }
      }

      k <- length(idx)
      lhs <- rbind(
        cbind(2.0 * risk_aversion * sigma, rep(1.0, k)),
        c(rep(1.0, k), 0.0)
      )
      rhs <- c(e, max_weight)
      sol <- try(qr.solve(lhs, rhs), silent = TRUE)
      if (!inherits(sol, "try-error")) {
        try_candidate(idx, as.numeric(sol[seq_len(k)]))
      }
    }
  }

  list(weights = best_w, objective = best_obj)
}

long_only_max_sharpe_opt <- function(edge,
                                     covmat,
                                     max_weight = 1.0,
                                     tol = 1e-10) {
  n <- length(edge)
  best_w <- numeric(n)
  best_sharpe <- -Inf

  sharpe_ratio <- function(w) {
    num <- sum(edge * w)
    den <- sqrt(as.numeric(t(w) %*% covmat %*% w))
    if (!is.finite(den) || den <= tol) {
      return(-Inf)
    }
    num / den
  }

  try_candidate <- function(idx, w_active) {
    if (any(!is.finite(w_active))) {
      return()
    }
    if (any(w_active < -tol)) {
      return()
    }

    s <- sum(w_active)
    if (s <= tol) {
      return()
    }

    w <- numeric(n)
    w[idx] <- max_weight * pmax(w_active, 0.0) / sum(pmax(w_active, 0.0))
    sr <- sharpe_ratio(w)

    if (is.finite(sr) && sr > best_sharpe + 1e-12) {
      best_w <<- w
      best_sharpe <<- sr
    }
  }

  for (active_size in seq_len(n)) {
    active_sets <- combn(seq_len(n), active_size, simplify = FALSE)

    for (idx in active_sets) {
      sigma <- covmat[idx, idx, drop = FALSE]
      e <- edge[idx]

      w_tan <- try(qr.solve(sigma, e), silent = TRUE)
      if (!inherits(w_tan, "try-error")) {
        try_candidate(idx, as.numeric(w_tan))
      }
    }
  }

  list(weights = best_w, sharpe = best_sharpe)
}

long_only_max_adjusted_sharpe_opt <- function(m,
                                              v,
                                              K,
                                              type,
                                              prices,
                                              edge,
                                              covmat,
                                              payoff_scenarios = NULL,
                                              max_weight = 1.0,
                                              tol = 1e-10) {
  n <- length(edge)

  softmax <- function(theta) {
    z <- theta - max(theta)
    ez <- exp(z)
    max_weight * ez / sum(ez)
  }

  adjusted_sharpe_for_weights <- function(w) {
    contracts <- w / prices
    stats <- if (is.null(payoff_scenarios)) {
      option_portfolio_stats(m, v, K, type, contracts)
    } else {
      option_portfolio_stats_from_scenarios(payoff_scenarios, contracts)
    }
    sharpe <- sum(edge * w) / stats["sd"]
    as.numeric(modified_sharpe_ratio(sharpe, stats["skew"], stats["ex.kurt"]))
  }

  objective <- function(theta) {
    val <- adjusted_sharpe_for_weights(softmax(theta))
    if (!is.finite(val)) {
      return(1e100)
    }
    -val
  }

  starts <- list(rep(0.0, n))
  for (i in seq_len(n)) {
    z <- rep(-8.0, n)
    z[i] <- 8.0
    starts[[length(starts) + 1]] <- z
  }

  sharpe_start <- long_only_max_sharpe_opt(edge, covmat, max_weight, tol)$weights
  if (sum(sharpe_start) > tol) {
    starts[[length(starts) + 1]] <- log(pmax(sharpe_start, tol))
  }

  best_w <- rep(max_weight / n, n)
  best_adj_sharpe <- adjusted_sharpe_for_weights(best_w)

  for (start in starts) {
    fit <- try(optim(start, objective, method = "BFGS", control = list(maxit = 1000)), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      w <- softmax(fit$par)
      val <- adjusted_sharpe_for_weights(w)
      if (is.finite(val) && val > best_adj_sharpe + 1e-12) {
        best_w <- w
        best_adj_sharpe <- val
      }
    }
  }

  list(weights = best_w, adjusted_sharpe = best_adj_sharpe)
}

signed_weight_contracts <- function(weights, bid_prices, ask_prices, budget) {
  exec_prices <- ifelse(weights < 0.0, bid_prices, ask_prices)
  budget * weights / exec_prices
}

long_short_portfolio_opt <- function(optimization_objective,
                                     m,
                                     v,
                                     K,
                                     type,
                                     bid_prices,
                                     ask_prices,
                                     budget,
                                     rf_growth,
                                     expected_payoff,
                                     cov_payoff,
                                     payoff_stats,
                                     payoff_scenarios,
                                     payoff_moment_scenarios,
                                     payoff_grid,
                                     tail_slope,
                                     risk_aversion,
                                     dependency_matrix = NULL,
                                     gamma = NA_real_,
                                     initial_contracts = rep(0.0, length(K)),
                                     max_gross_weight = 1.0,
                                     force_full_investment = TRUE,
                                     tol = 1e-10) {
  n <- length(K)
  theta_len <- if (force_full_investment) 2L * n else 2L * n + 1L

  softmax <- function(theta) {
    z <- theta - max(theta)
    ez <- exp(z)
    ez / sum(ez)
  }

  weights_from_theta <- function(theta) {
    sleeve_weights <- max_gross_weight * softmax(theta)
    if (!force_full_investment) {
      sleeve_weights <- sleeve_weights[-length(sleeve_weights)]
    }
    sleeve_weights[seq_len(n)] - sleeve_weights[n + seq_len(n)]
  }

  objective_value <- function(weights) {
    contracts <- signed_weight_contracts(weights, bid_prices, ask_prices, budget)
    portfolio_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = contracts,
      initial_contracts = initial_contracts,
      bid_prices = bid_prices,
      ask_prices = ask_prices,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      min_terminal_wealth = -Inf,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix
    )
    portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = portfolio_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = ask_prices,
      budget = budget,
      rf_growth = rf_growth,
      contracts = contracts,
      utility_gamma = gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_prices,
      ask_prices = ask_prices,
      payoff_scenarios = payoff_scenarios
    )$value
  }

  objective <- function(theta) {
    val <- objective_value(weights_from_theta(theta))
    if (!isTRUE(is.finite(val))) {
      return(1e100)
    }
    -val
  }

  long_edge <- expected_payoff / ask_prices - rf_growth
  short_edge <- rf_growth - expected_payoff / bid_prices
  starts <- list(rep(0.0, theta_len))
  long_idx <- which.max(long_edge)
  short_idx <- which.max(short_edge)
  if (is.finite(long_edge[long_idx])) {
    long_start <- rep(-6.0, theta_len)
    long_start[long_idx] <- 6.0
    starts[[length(starts) + 1L]] <- long_start
  }
  if (is.finite(short_edge[short_idx])) {
    short_start <- rep(-6.0, theta_len)
    short_start[n + short_idx] <- 6.0
    starts[[length(starts) + 1L]] <- short_start
  }
  if (is.finite(long_edge[long_idx]) && is.finite(short_edge[short_idx])) {
    pair_start <- rep(-6.0, theta_len)
    pair_start[long_idx] <- 6.0
    pair_start[n + short_idx] <- 6.0
    starts[[length(starts) + 1L]] <- pair_start
  }

  best_theta <- starts[[1L]]
  best_weights <- weights_from_theta(best_theta)
  best_objective <- objective_value(best_weights)
  if (!isTRUE(is.finite(best_objective))) {
    best_objective <- -Inf
  }
  best_optimizer <- "initial"
  best_status <- "initial"
  best_status_code <- NA_integer_
  best_evals <- 0L

  for (start in starts) {
    fit <- try(optim(
      start,
      objective,
      method = "BFGS",
      control = list(maxit = 80)
    ), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      weights <- weights_from_theta(fit$par)
      val <- objective_value(weights)
      if (isTRUE(is.finite(val)) && val > best_objective + 1e-12) {
        best_theta <- fit$par
        best_weights <- weights
        best_objective <- val
        best_optimizer <- "optim_bfgs"
        best_status <- decode_optim_status(fit$convergence)
        best_status_code <- as.integer(fit$convergence)
        best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
      }
    }
  }

  contracts <- signed_weight_contracts(best_weights, bid_prices, ask_prices, budget)
  list(
    weights = best_weights,
    contracts = contracts,
    trade_contracts = rebalance_trade_contracts(contracts, initial_contracts),
    objective = best_objective,
    optimizer_diagnostics = list(
      optimizer = best_optimizer,
      optimizer_status = best_status,
      optimizer_status_code = best_status_code,
      optimizer_evals = best_evals,
      optimizer_fallback_used = FALSE
    )
  )
}

option_portfolio_raw_moment <- function(m, v, K, type, weights, n) {
  stopifnot(length(K) == length(type))
  stopifnot(length(K) == length(weights))

  stops <- sort(unique(c(0, K)))
  total <- 0.0

  for (j in seq_along(stops)) {
    lower <- stops[j]
    if (j < length(stops)) {
      upper <- stops[j + 1]
    } else {
      upper <- Inf
    }

    call_active <- type == "call" & K <= lower
    put_active <- type == "put" & K >= upper

    a <- sum(weights[call_active]) - sum(weights[put_active])
    b <- -sum(weights[call_active] * K[call_active]) +
      sum(weights[put_active] * K[put_active])

    interval_value <- 0.0
    for (ell in 0:n) {
      interval_value <- interval_value +
        choose(n, ell) *
        a^ell *
        b^(n - ell) *
        ln_interval_moment(m, v, lower, upper, ell)
    }

    total <- total + interval_value
  }

  total
}

option_portfolio_stats <- function(m, v, K, type, weights) {
  r1 <- option_portfolio_raw_moment(m, v, K, type, weights, 1)
  r2 <- option_portfolio_raw_moment(m, v, K, type, weights, 2)
  r3 <- option_portfolio_raw_moment(m, v, K, type, weights, 3)
  r4 <- option_portfolio_raw_moment(m, v, K, type, weights, 4)

  mu2 <- r2 - r1^2
  mu3 <- r3 - 3 * r1 * r2 + 2 * r1^3
  mu4 <- r4 - 4 * r1 * r3 + 6 * r1^2 * r2 - 3 * r1^4

  c(
    mean = r1,
    sd = sqrt(mu2),
    skew = mu3 / mu2^(3 / 2),
    ex.kurt = mu4 / mu2^2 - 3
  )
}

option_payoff_stats_table <- function(m, v, K, type) {
  n <- length(K)
  means <- numeric(n)
  sds <- numeric(n)
  skews <- numeric(n)
  exkurts <- numeric(n)

  for (i in seq_len(n)) {
    stats <- option_portfolio_stats(m, v, K[i], type[i], 1.0)
    means[i] <- stats["mean"]
    sds[i] <- stats["sd"]
    skews[i] <- stats["skew"]
    exkurts[i] <- stats["ex.kurt"]
  }

  data.frame(
    type = type,
    strike = K,
    mean = means,
    sd = sds,
    skew = skews,
    ex.kurt = exkurts
  )
}

option_payoff_cov_mat <- function(m, v, K, type) {
  n <- length(K)
  out <- matrix(0.0, n, n)
  second <- numeric(n)

  for (i in seq_len(n)) {
    second[i] <- option_portfolio_raw_moment(m, v, K[i], type[i], 1.0, 2)
  }

  means <- option_payoff_stats_table(m, v, K, type)$mean

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      e_sum_sq <- option_portfolio_raw_moment(
        m,
        v,
        c(K[i], K[j]),
        c(type[i], type[j]),
        c(1.0, 1.0),
        2
      )
      eij <- 0.5 * (e_sum_sq - second[i] - second[j])
      out[i, j] <- eij - means[i] * means[j]
    }
  }

  labels <- paste(type, K, sep = "=")
  rownames(out) <- labels
  colnames(out) <- labels
  out
}

sample_moments <- function(x) {
  mean_x <- mean(x)
  centered <- x - mean_x
  mu2 <- mean(centered^2)

  if (!is.finite(mu2) || mu2 <= 0.0) {
    return(c(mean = mean_x, sd = 0.0, skew = NaN, ex.kurt = NaN))
  }

  mu3 <- mean(centered^3)
  mu4 <- mean(centered^4)
  c(
    mean = mean_x,
    sd = sqrt(mu2),
    skew = mu3 / mu2^(3 / 2),
    ex.kurt = mu4 / mu2^2 - 3.0
  )
}

option_payoff_stats_table_from_scenarios <- function(payoff_scenarios, K, type) {
  n <- length(K)
  means <- numeric(n)
  sds <- numeric(n)
  skews <- numeric(n)
  exkurts <- numeric(n)

  for (i in seq_len(n)) {
    stats <- sample_moments(payoff_scenarios[, i])
    means[i] <- stats["mean"]
    sds[i] <- stats["sd"]
    skews[i] <- stats["skew"]
    exkurts[i] <- stats["ex.kurt"]
  }

  data.frame(
    type = type,
    strike = K,
    mean = means,
    sd = sds,
    skew = skews,
    ex.kurt = exkurts
  )
}

option_payoff_cov_mat_from_scenarios <- function(payoff_scenarios, K, type) {
  out <- cov(payoff_scenarios) * ((nrow(payoff_scenarios) - 1) / nrow(payoff_scenarios))
  labels <- paste(type, K, sep = "=")
  rownames(out) <- labels
  colnames(out) <- labels
  out
}

option_portfolio_stats_from_scenarios <- function(payoff_scenarios, weights) {
  sample_moments(as.numeric(payoff_scenarios %*% weights))
}

option_payoff_matrix <- function(S, K, type) {
  out <- matrix(0.0, length(S), length(K))

  for (i in seq_along(K)) {
    if (type[i] == "call") {
      out[, i] <- pmax(S - K[i], 0.0)
    } else if (type[i] == "put") {
      out[, i] <- pmax(K[i] - S, 0.0)
    } else {
      stop("type must be 'call' or 'put'")
    }
  }

  out
}

multivariate_lognormal_terminal_prices <- function(S0,
                                                   mu,
                                                   cov_matrix,
                                                   T,
                                                   n_scenarios = 1001,
                                                   seed = NA_integer_) {
  S0 <- as.numeric(S0)
  mu <- as.numeric(mu)
  stopifnot(length(S0) == length(mu))
  stopifnot(nrow(cov_matrix) == length(S0))
  stopifnot(ncol(cov_matrix) == length(S0))
  if (!is.na(seed)) {
    set.seed(as.integer(seed))
  }
  n_assets <- length(S0)
  chol_cov <- chol(cov_matrix * T)
  z <- matrix(rnorm(n_scenarios * n_assets), nrow = n_scenarios, ncol = n_assets)
  shocks <- z %*% chol_cov
  variance <- diag(cov_matrix)
  drift <- matrix(
    log(S0) + (mu - 0.5 * variance) * T,
    nrow = n_scenarios,
    ncol = n_assets,
    byrow = TRUE
  )
  exp(drift + shocks)
}

multistock_payoff_matrix <- function(terminal_prices,
                                     underlying,
                                     K,
                                     type,
                                     underlying_names) {
  out <- matrix(0.0, nrow = nrow(terminal_prices), ncol = length(K))
  for (j in seq_along(K)) {
    asset_idx <- match(underlying[j], underlying_names)
    if (is.na(asset_idx)) {
      stop("Unknown underlying: ", underlying[j])
    }
    S <- terminal_prices[, asset_idx]
    if (K[j] <= 0.0 && type[j] == "call") {
      out[, j] <- S
    } else if (type[j] == "call") {
      out[, j] <- pmax(S - K[j], 0.0)
    } else if (type[j] == "put") {
      out[, j] <- pmax(K[j] - S, 0.0)
    } else {
      stop("type must be 'call' or 'put'")
    }
  }
  out
}

multistock_price_vec <- function(S0,
                                 underlying,
                                 K,
                                 type,
                                 r,
                                 q,
                                 sigma,
                                 T,
                                 underlying_names) {
  out <- numeric(length(K))
  for (j in seq_along(K)) {
    asset_idx <- match(underlying[j], underlying_names)
    if (is.na(asset_idx)) {
      stop("Unknown underlying: ", underlying[j])
    }
    if (K[j] <= 0.0 && type[j] == "call") {
      out[j] <- S0[asset_idx]
    } else {
      out[j] <- bs_option_price_vec(S0[asset_idx], K[j], r, q[asset_idx], sigma[j], T, type[j])
    }
  }
  out
}

multistock_greeks_table <- function(S0,
                                    underlying,
                                    K,
                                    type,
                                    r,
                                    q,
                                    sigma,
                                    T,
                                    underlying_names) {
  out <- t(vapply(seq_along(K), function(j) {
    asset_idx <- match(underlying[j], underlying_names)
    if (is.na(asset_idx)) {
      stop("Unknown underlying: ", underlying[j])
    }
    if (K[j] <= 0.0 && type[j] == "call") {
      c(delta = exp(-q[asset_idx] * T), gamma = 0.0, vega = 0.0)
    } else {
      bs_option_greeks(S0[asset_idx], K[j], r, q[asset_idx], sigma[j], T, type[j])
    }
  }, numeric(3)))
  out <- as.data.frame(out)
  rownames(out) <- NULL
  out
}

multistock_horizon_inputs <- function(terminal_prices,
                                      cvar_prices,
                                      underlying,
                                      K,
                                      type,
                                      underlying_names) {
  payoff_scenarios <- multistock_payoff_matrix(
    terminal_prices = terminal_prices,
    underlying = underlying,
    K = K,
    type = type,
    underlying_names = underlying_names
  )
  cvar_payoff_scenarios <- multistock_payoff_matrix(
    terminal_prices = cvar_prices,
    underlying = underlying,
    K = K,
    type = type,
    underlying_names = underlying_names
  )
  payoff_stats <- option_payoff_stats_table_from_scenarios(payoff_scenarios, K, type)
  list(
    payoff_scenarios = payoff_scenarios,
    cvar_payoff_scenarios = cvar_payoff_scenarios,
    payoff_stats = payoff_stats,
    expected_payoff = payoff_stats$mean,
    cov_payoff = option_payoff_cov_mat_from_scenarios(payoff_scenarios, K, type)
  )
}

multitenor_option_value_matrix <- function(S,
                                           K,
                                           type,
                                           expiry,
                                           horizon,
                                           r,
                                           q,
                                           sigma) {
  if (length(sigma) == 1L) {
    sigma <- rep(sigma, length(K))
  }
  stopifnot(length(K) == length(type))
  stopifnot(length(K) == length(expiry))
  stopifnot(length(K) == length(sigma))

  out <- matrix(0.0, nrow = length(S), ncol = length(K))
  for (j in seq_along(K)) {
    if (K[j] <= 0.0 && type[j] == "call") {
      out[, j] <- S
    } else if (expiry[j] <= horizon + 1e-12) {
      out[, j] <- option_payoff_matrix(S, K[j], type[j])[, 1L]
    } else {
      remaining_T <- expiry[j] - horizon
      out[, j] <- vapply(S, function(s) {
        bs_option_price_vec(s, K[j], r, q, sigma[j], remaining_T, type[j])
      }, numeric(1))
    }
  }
  out
}

multitenor_price_vec <- function(S0,
                                 K,
                                 type,
                                 expiry,
                                 r,
                                 q,
                                 sigma) {
  if (length(sigma) == 1L) {
    sigma <- rep(sigma, length(K))
  }
  out <- numeric(length(K))
  for (j in seq_along(K)) {
    if (K[j] <= 0.0 && type[j] == "call") {
      out[j] <- S0
    } else {
      out[j] <- bs_option_price_vec(S0, K[j], r, q, sigma[j], expiry[j], type[j])
    }
  }
  out
}

multitenor_greeks_table <- function(S0,
                                    K,
                                    type,
                                    expiry,
                                    r,
                                    q,
                                    sigma) {
  if (length(sigma) == 1L) {
    sigma <- rep(sigma, length(K))
  }
  out <- t(vapply(seq_along(K), function(j) {
    bs_option_greeks(S0, K[j], r, q, sigma[j], expiry[j], type[j])
  }, numeric(3)))
  out <- as.data.frame(out)
  rownames(out) <- NULL
  out
}

multitenor_horizon_inputs <- function(horizon_prices,
                                      cvar_prices,
                                      state_grid,
                                      K,
                                      type,
                                      expiry,
                                      horizon,
                                      r,
                                      q,
                                      sigma) {
  value_scenarios <- multitenor_option_value_matrix(
    S = horizon_prices,
    K = K,
    type = type,
    expiry = expiry,
    horizon = horizon,
    r = r,
    q = q,
    sigma = sigma
  )
  cvar_value_scenarios <- multitenor_option_value_matrix(
    S = cvar_prices,
    K = K,
    type = type,
    expiry = expiry,
    horizon = horizon,
    r = r,
    q = q,
    sigma = sigma
  )
  value_grid <- multitenor_option_value_matrix(
    S = state_grid,
    K = K,
    type = type,
    expiry = expiry,
    horizon = horizon,
    r = r,
    q = q,
    sigma = sigma
  )
  value_stats <- option_payoff_stats_table_from_scenarios(value_scenarios, K, type)
  list(
    value_scenarios = value_scenarios,
    cvar_value_scenarios = cvar_value_scenarios,
    value_grid = value_grid,
    value_stats = value_stats,
    expected_value = value_stats$mean,
    cov_value = option_payoff_cov_mat_from_scenarios(value_scenarios, K, type)
  )
}

option_tail_slope <- function(K, type) {
  ifelse(type == "call", 1.0, 0.0)
}

option_lower_tail_slope <- function(K, type) {
  ifelse(type == "put", -1.0, 0.0)
}

terminal_state_grid <- function(K, terminal_floor = NA_real_, lower_unbounded = FALSE) {
  if (lower_unbounded) {
    return(sort(unique(K)))
  }
  if (!is.na(terminal_floor)) {
    return(sort(unique(c(terminal_floor, K[K >= terminal_floor]))))
  }
  sort(unique(c(0.0, K)))
}

terminal_distribution_inputs <- function(terminal_model,
                                         S0,
                                         mu,
                                         sigma,
                                         T,
                                         K,
                                         type,
                                         terminal_floor = NA_real_,
                                         n_scenarios = 401,
                                         cvar_n_scenarios = 401,
                                         mixture_weights = NULL,
                                         mixture_mu = NULL,
                                         mixture_sigma = NULL,
                                         hyperbolic_alpha = 10.0,
                                         hyperbolic_beta = -2.0,
                                         gh_lambda = 1.0,
                                         gh_alpha = 10.0,
                                         gh_beta = -2.0,
                                         nig_alpha = 10.0,
                                         nig_beta = -2.0,
                                         vg_shape = 1.0,
                                         vg_beta = -0.2) {
  section_timings <- data.frame(section = character(), elapsed_seconds = numeric())
  time_terminal_section <- function(section, expr) {
    section_start <- proc.time()
    value <- force(expr)
    section_elapsed <- proc.time() - section_start
    section_timings <<- rbind(
      section_timings,
      data.frame(section = section, elapsed_seconds = as.numeric(section_elapsed["elapsed"]))
    )
    value
  }

  par <- if (terminal_model == "lognormal") {
    ln_params_from_gbm(S0, mu, sigma, T)
  } else if (terminal_model == "logistic") {
    logistic_par <- logistic_logreturn_params(S0, mu, sigma, T)
    list(m = logistic_par$location, v = logistic_par$v)
  } else if (terminal_model == "hyperbolic_secant") {
    hs_par <- hyperbolic_secant_logreturn_params(S0, mu, sigma, T)
    list(m = hs_par$location, v = hs_par$v)
  } else if (terminal_model == "symmetric_hyperbolic") {
    hyperbolic_par <- symmetric_hyperbolic_logreturn_params(S0, mu, sigma, T, alpha = hyperbolic_alpha)
    list(m = hyperbolic_par$location, v = hyperbolic_par$v)
  } else if (terminal_model == "hyperbolic") {
    hyperbolic_par <- hyperbolic_logreturn_params(S0, mu, sigma, T, alpha = hyperbolic_alpha, beta = hyperbolic_beta)
    list(m = hyperbolic_par$location, v = hyperbolic_par$v)
  } else if (terminal_model == "generalized_hyperbolic") {
    gh_par <- gh_logreturn_params(S0, mu, sigma, T, lambda = gh_lambda, alpha = gh_alpha, beta = gh_beta)
    list(m = gh_par$location, v = gh_par$v)
  } else if (terminal_model == "normal_inverse_gaussian") {
    nig_par <- gh_logreturn_params(S0, mu, sigma, T, lambda = -0.5, alpha = nig_alpha, beta = nig_beta,
                                   label = "normal_inverse_gaussian")
    list(m = nig_par$location, v = nig_par$v)
  } else if (terminal_model == "variance_gamma") {
    vg_par <- vg_logreturn_params(S0, mu, sigma, T, shape = vg_shape, beta = vg_beta)
    list(m = vg_par$location, v = vg_par$v)
  } else {
    list(m = NA_real_, v = NA_real_)
  }
  lower_tail_slope <- if (terminal_model == "normal" && is.na(terminal_floor)) option_lower_tail_slope(K, type) else NULL
  scenario_prices <- time_terminal_section("scenario_quantiles", terminal_price_quantiles(
    model = terminal_model,
    S0 = S0,
    mu = mu,
    sigma = sigma,
    T = T,
    n = n_scenarios,
    terminal_floor = terminal_floor,
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
  payoff_scenarios <- time_terminal_section("payoff_scenarios", option_payoff_matrix(scenario_prices, K, type))

  payoff_moments <- time_terminal_section("payoff_moments", {
    if (terminal_model == "lognormal" && is.na(terminal_floor)) {
      payoff_stats <- option_payoff_stats_table(par$m, par$v, K, type)
      expected_payoff <- payoff_stats$mean
      cov_payoff <- option_payoff_cov_mat(par$m, par$v, K, type)
    } else {
      payoff_stats <- option_payoff_stats_table_from_scenarios(payoff_scenarios, K, type)
      expected_payoff <- payoff_stats$mean
      cov_payoff <- option_payoff_cov_mat_from_scenarios(payoff_scenarios, K, type)
    }
    list(payoff_stats = payoff_stats, expected_payoff = expected_payoff, cov_payoff = cov_payoff)
  })
  payoff_stats <- payoff_moments$payoff_stats
  expected_payoff <- payoff_moments$expected_payoff
  cov_payoff <- payoff_moments$cov_payoff

  grid_inputs <- time_terminal_section("terminal_constraint_grid", {
    state_grid <- terminal_state_grid(K, terminal_floor, lower_unbounded = !is.null(lower_tail_slope))
    payoff_grid <- option_payoff_matrix(state_grid, K, type)
    list(state_grid = state_grid, payoff_grid = payoff_grid)
  })
  state_grid <- grid_inputs$state_grid
  payoff_grid <- grid_inputs$payoff_grid
  cvar_scenarios <- time_terminal_section("cvar_quantiles", terminal_price_quantiles(
    model = terminal_model,
    S0 = S0,
    mu = mu,
    sigma = sigma,
    T = T,
    n = cvar_n_scenarios,
    terminal_floor = terminal_floor,
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
  cvar_payoff_scenarios <- time_terminal_section("cvar_payoff_scenarios", option_payoff_matrix(cvar_scenarios, K, type))

  list(
    m = par$m,
    v = par$v,
    scenario_prices = scenario_prices,
    payoff_scenarios = payoff_scenarios,
    payoff_stats = payoff_stats,
    expected_payoff = expected_payoff,
    cov_payoff = cov_payoff,
    state_grid = state_grid,
    payoff_grid = payoff_grid,
    tail_slope = option_tail_slope(K, type),
    lower_tail_slope = lower_tail_slope,
    cvar_scenarios = cvar_scenarios,
    cvar_payoff_scenarios = cvar_payoff_scenarios,
    section_timings = section_timings
  )
}

cvar_loss_values <- function(terminal_wealth, budget, tail_probs) {
  loss <- budget - terminal_wealth
  sorted_loss <- sort(loss, decreasing = TRUE)
  n <- length(sorted_loss)
  out <- numeric(length(tail_probs))

  for (i in seq_along(tail_probs)) {
    k <- max(1, ceiling(tail_probs[i] * n))
    out[i] <- mean(sorted_loss[seq_len(k)])
  }

  out
}

cvar_constraint_penalty <- function(terminal_wealth,
                                    budget,
                                    cvar_constraints,
                                    penalty_scale = 1e6) {
  if (is.null(cvar_constraints) || nrow(cvar_constraints) == 0) {
    return(0.0)
  }

  es <- cvar_loss_values(terminal_wealth, budget, cvar_constraints$tail_prob)
  excess <- pmax(es - cvar_constraints$max_loss, 0.0)
  penalty_scale * sum(excess^2) / max(1.0, budget^2)
}

cvar_constraint_excess <- function(terminal_wealth, budget, cvar_constraints) {
  if (is.null(cvar_constraints) || nrow(cvar_constraints) == 0) {
    return(numeric(0))
  }

  cvar_loss_values(terminal_wealth, budget, cvar_constraints$tail_prob) -
    cvar_constraints$max_loss
}

cvar_constraint_jacobian <- function(contracts,
                                     budget,
                                     rf_growth,
                                     cvar_constraints,
                                     cvar_payoff_scenarios,
                                     initial_contracts,
                                     bid_prices,
                                     ask_prices) {
  if (is.null(cvar_constraints) || nrow(cvar_constraints) == 0 ||
      is.null(cvar_payoff_scenarios) || nrow(cvar_payoff_scenarios) == 0) {
    return(matrix(0.0, 0, length(contracts)))
  }

  trade_cost <- rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
  cash_growth <- (budget - trade_cost) * rf_growth
  terminal_wealth <- cash_growth + as.numeric(cvar_payoff_scenarios %*% contracts)
  loss_order <- order(budget - terminal_wealth, decreasing = TRUE)
  trade_grad <- position_exec_prices(
    rebalance_trade_contracts(contracts, initial_contracts),
    bid_prices,
    ask_prices
  )
  jac <- matrix(0.0, nrow(cvar_constraints), length(contracts))

  for (i in seq_len(nrow(cvar_constraints))) {
    k <- max(1, ceiling(cvar_constraints$tail_prob[i] * length(loss_order)))
    rows <- loss_order[seq_len(k)]
    jac[i, ] <- colMeans(
      rf_growth * matrix(trade_grad, k, length(contracts), byrow = TRUE) -
        cvar_payoff_scenarios[rows, , drop = FALSE]
    )
  }

  jac
}

es_loss_col_names <- function(tail_probs) {
  pct <- 100 * tail_probs
  pct_label <- ifelse(abs(pct - round(pct)) < 1e-8,
                      as.character(as.integer(round(pct))),
                      gsub("\\.", "p", format(pct, trim = TRUE, scientific = FALSE)))
  paste0("es_loss_", pct_label, "pct")
}

safe_file_component <- function(x) {
  x <- ifelse(is.na(x), "na", as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

option_run_table_csv_file <- function(run, dir) {
  parts <- c(run$constraint_mode, run$objective)
  if (!is.na(run$utility_gamma)) {
    gamma <- format(run$utility_gamma, trim = TRUE, scientific = FALSE)
    parts <- c(parts, paste0("gamma_", gamma))
  }
  file.path(dir, paste0(paste(safe_file_component(parts), collapse = "_"), ".csv"))
}

write_option_portfolio_tables_csv <- function(runs, dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  files <- character(length(runs))
  for (i in seq_along(runs)) {
    file <- option_run_table_csv_file(runs[[i]], dir)
    table <- runs[[i]]$table_raw
    table$constraint_mode <- runs[[i]]$constraint_mode
    table$objective <- runs[[i]]$objective
    table$utility_gamma <- runs[[i]]$utility_gamma
    write.csv(table, file, row.names = FALSE)
    files[i] <- file
  }
  files
}

combined_option_portfolio_table <- function(runs) {
  tables <- vector("list", length(runs))
  for (i in seq_along(runs)) {
    table <- runs[[i]]$table_raw
    table$run_id <- i
    table$constraint_mode <- runs[[i]]$constraint_mode
    table$objective <- runs[[i]]$objective
    table$utility_gamma <- runs[[i]]$utility_gamma
    tables[[i]] <- table
  }
  do.call(rbind, tables)
}

write_combined_option_portfolio_table_csv <- function(runs, file) {
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  table <- combined_option_portfolio_table(runs)
  write.csv(table, file, row.names = FALSE)
  invisible(file)
}

normalize_portfolio_delta_bounds <- function(bounds) {
  if (is.null(bounds) || length(bounds) == 0) {
    return(c(NA_real_, NA_real_))
  }
  if (length(bounds) != 2) {
    stop("portfolio_delta_bounds must have length 2: c(lower, upper)")
  }
  bounds <- as.numeric(bounds)
  if (is.na(bounds[1])) bounds[1] <- -Inf
  if (is.na(bounds[2])) bounds[2] <- Inf
  if (bounds[1] > bounds[2]) {
    stop("portfolio_delta_bounds lower bound exceeds upper bound")
  }
  bounds
}

portfolio_delta_constraint_rows <- function(delta, bounds) {
  bounds <- normalize_portfolio_delta_bounds(bounds)
  if (is.null(delta) && any(is.finite(bounds))) {
    stop("portfolio_delta_bounds require report_greeks = TRUE so deltas are available")
  }
  ui <- NULL
  ci <- numeric(0)
  if (is.finite(bounds[1])) {
    ui <- rbind(ui, delta)
    ci <- c(ci, bounds[1])
  }
  if (is.finite(bounds[2])) {
    ui <- rbind(ui, -delta)
    ci <- c(ci, -bounds[2])
  }
  list(ui = ui, ci = ci, bounds = bounds)
}

portfolio_delta_bound_violation <- function(portfolio_delta, bounds) {
  bounds <- normalize_portfolio_delta_bounds(bounds)
  if (!is.finite(portfolio_delta)) {
    return(NA_real_)
  }
  max(bounds[1] - portfolio_delta, portfolio_delta - bounds[2], 0.0)
}

scale_contracts_to_terminal_constraints <- function(contracts,
                                                    prices,
                                                    budget,
                                                    rf_growth,
                                                    payoff_grid,
                                                    lower_tail_slope = NULL,
                                                    min_terminal_wealth = -Inf,
                                                    cvar_constraints = NULL,
                                                    cvar_payoff_scenarios = NULL,
                                                    delta = NULL,
                                                    portfolio_delta_bounds = c(NA_real_, NA_real_),
                                                    initial_contracts = rep(0.0, length(contracts)),
                                                    bid_prices = prices,
                                                    ask_prices = prices,
                                                    tol = 1e-8,
                                                    max_iter = 80) {
  feasible <- function(lambda) {
    scaled_contracts <- initial_contracts + lambda * (contracts - initial_contracts)
    cash <- budget - rebalance_trade_cost(
      scaled_contracts,
      initial_contracts,
      bid_prices,
      ask_prices
    )
    terminal_grid <- cash * rf_growth + as.numeric(payoff_grid %*% scaled_contracts)

    if (min(terminal_grid) < min_terminal_wealth - tol) {
      return(FALSE)
    }

    if (!is.null(lower_tail_slope) && sum(lower_tail_slope * scaled_contracts) > tol) {
      return(FALSE)
    }

    if (!is.null(delta)) {
      if (portfolio_delta_bound_violation(sum(delta * scaled_contracts), portfolio_delta_bounds) > tol) {
        return(FALSE)
      }
    }

    if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) {
      terminal_cvar <- cash * rf_growth +
        as.numeric(cvar_payoff_scenarios %*% scaled_contracts)
      if (any(cvar_constraint_excess(terminal_cvar, budget, cvar_constraints) > tol)) {
        return(FALSE)
      }
    }

    TRUE
  }

  if (feasible(1.0)) {
    return(contracts)
  }

  lo <- 0.0
  hi <- 1.0
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (feasible(mid)) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  initial_contracts + lo * (contracts - initial_contracts)
}

expected_utility_for_contracts <- function(m,
                                           v,
                                           K,
                                           type,
                                           prices,
                                           budget,
                                           rf_growth,
                                           contracts,
                                           gamma,
                                           initial_contracts = rep(0.0, length(contracts)),
                                           bid_prices = prices,
                                           ask_prices = prices,
                                           payoff_scenarios = NULL,
                                           n_grid = 401) {
  if (is.null(payoff_scenarios)) {
    z <- qnorm(seq(0.001, 0.999, length.out = n_grid))
    S <- exp(m + sqrt(v) * z)
    payoff_scenarios <- option_payoff_matrix(S, K, type)
  }
  cash <- budget - rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
  terminal_wealth <- cash * rf_growth + as.numeric(payoff_scenarios %*% contracts)
  expected_utility_value(terminal_wealth, gamma)
}

evaluate_option_portfolio <- function(m,
                                      v,
                                      K,
                                      type,
                                      contracts,
                                      initial_contracts,
                                      bid_prices,
                                      ask_prices,
                                      budget,
                                      rf_growth,
                                      expected_payoff,
                                      cov_payoff,
                                      payoff_stats,
                                      payoff_grid,
                                      payoff_scenarios = NULL,
                                      tail_slope,
                                      cvar_constraints = NULL,
                                      cvar_payoff_scenarios = NULL,
                                      min_terminal_wealth = -Inf,
                                      risk_aversion = 0.0,
                                      dependency_matrix = NULL,
                                      greeks = NULL,
                                      S0 = NA_real_,
                                      portfolio_delta_bounds = c(NA_real_, NA_real_),
                                      mtm_var_inputs = NULL,
                                      mtm_var_constraints = NULL,
                                      mtm_es_constraints = NULL,
                                      delta_hedge_inputs = NULL) {
  trade_contracts <- rebalance_trade_contracts(contracts, initial_contracts)
  trade_cost <- executable_trade_cost(trade_contracts, bid_prices, ask_prices)
  cash_after_trade <- budget - trade_cost
  mean_wealth <- cash_after_trade * rf_growth + sum(expected_payoff * contracts)
  variance_wealth <- as.numeric(t(contracts) %*% cov_payoff %*% contracts)
  sd_wealth <- sqrt(max(variance_wealth, 0.0))
  edge_wealth <- mean_wealth - budget * rf_growth
  sharpe <- edge_wealth / sd_wealth

  portfolio_cov <- as.numeric(cov_payoff %*% contracts)
  port_corr <- portfolio_cov / (payoff_stats$sd * sd_wealth)
  payoff_moments <- if (is.null(payoff_scenarios)) {
    option_portfolio_stats(m, v, K, type, contracts)
  } else {
    option_portfolio_stats_from_scenarios(payoff_scenarios, contracts)
  }
  adj_sharpe <- modified_sharpe_ratio(
    sharpe,
    payoff_moments["skew"],
    payoff_moments["ex.kurt"]
  )

  min_wealth <- min(cash_after_trade * rf_growth + as.numeric(payoff_grid %*% contracts))
  min_wealth_violation <- max(min_terminal_wealth - min_wealth, 0.0)
  tail_slope_value <- sum(tail_slope * contracts)

  es_loss <- numeric(0)
  es_names <- character(0)
  es_violation <- numeric(0)
  es_violation_names <- character(0)
  max_es_violation <- 0.0
  if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) {
    cvar_terminal_wealth <- cash_after_trade * rf_growth +
      as.numeric(cvar_payoff_scenarios %*% contracts)
    es_loss <- cvar_loss_values(cvar_terminal_wealth, budget, cvar_constraints$tail_prob)
    es_names <- es_loss_col_names(cvar_constraints$tail_prob)
    es_violation <- pmax(es_loss - cvar_constraints$max_loss, 0.0)
    es_violation_names <- sub("^es_loss_", "es_violation_", es_names)
    max_es_violation <- max(es_violation)
  }

  position_price <- position_exec_prices(contracts, bid_prices, ask_prices)
  trade_price <- position_exec_prices(trade_contracts, bid_prices, ask_prices)
  position_value <- position_price * contracts
  position_weight <- position_value / budget
  gross_position_weight <- sum(abs(position_weight))
  trade_cost_by_instrument <- trade_price * trade_contracts
  trade_weight <- trade_cost_by_instrument / budget

  payoff_per_dollar <- expected_payoff / trade_price
  edge <- payoff_per_dollar - rf_growth
  sd_per_dollar <- payoff_stats$sd / trade_price
  instrument_sharpe <- edge / sd_per_dollar
  instrument_adj_sharpe <- modified_sharpe_ratio(
    instrument_sharpe,
    payoff_stats$skew,
    payoff_stats$ex.kurt
  )
  dependency_objective <- NA_real_
  dependency_penalty <- NA_real_
  if (!is.null(dependency_matrix)) {
    dependency_weights <- position_exec_prices(contracts, bid_prices, ask_prices) * contracts / budget
    dependency_penalty <- as.numeric(t(dependency_weights) %*% dependency_matrix %*% dependency_weights)
    dependency_objective <- edge_wealth / budget - risk_aversion * dependency_penalty
  }
  portfolio_delta <- NA_real_
  portfolio_gamma <- NA_real_
  portfolio_vega <- NA_real_
  delta_dollars <- NA_real_
  gamma_dollars_1pct <- NA_real_
  vega_1pct <- NA_real_
  delta_bound_violation <- NA_real_
  if (!is.null(greeks)) {
    portfolio_delta <- sum(contracts * greeks$delta)
    portfolio_gamma <- sum(contracts * greeks$gamma)
    portfolio_vega <- sum(contracts * greeks$vega)
    delta_bound_violation <- portfolio_delta_bound_violation(
      portfolio_delta,
      portfolio_delta_bounds
    )
    if (is.finite(S0)) {
      delta_dollars <- S0 * portfolio_delta
      gamma_dollars_1pct <- 0.5 * portfolio_gamma * (0.01 * S0)^2
    }
    vega_1pct <- 0.01 * portfolio_vega
  }
  mtm_var <- list(
    current_mtm_wealth = NA_real_,
    mean_mtm_pnl = NA_real_,
    sd_mtm_pnl = NA_real_
  )
  if (!is.null(mtm_var_inputs)) {
    if (identical(mtm_var_inputs$mode, "multistock")) {
      mtm_var <- multistock_portfolio_mtm_var_stats(
        contracts = contracts,
        cash = cash_after_trade,
        current_mid_price = mtm_var_inputs$current_mid_price,
        S0 = mtm_var_inputs$S0,
        underlying = mtm_var_inputs$underlying,
        K = K,
        r = mtm_var_inputs$r,
        q = mtm_var_inputs$q,
        mid_vol = mtm_var_inputs$mid_vol,
        T = mtm_var_inputs$T,
        type = type,
        var_mu = mtm_var_inputs$var_mu,
        var_cov_matrix = mtm_var_inputs$var_cov_matrix,
        underlying_names = mtm_var_inputs$underlying_names,
        horizon_days = mtm_var_inputs$horizon_days,
        trading_days = mtm_var_inputs$trading_days,
        n_scenarios = mtm_var_inputs$n_scenarios,
        return_model = mtm_var_inputs$return_model,
        terminal_floor = mtm_var_inputs$terminal_floor,
        conf_levels = mtm_var_inputs$conf_levels,
        seed = mtm_var_inputs$seed
      )
    } else {
      mtm_var <- portfolio_mtm_var_stats(
        contracts = contracts,
        cash = cash_after_trade,
        current_mid_price = mtm_var_inputs$current_mid_price,
        S0 = mtm_var_inputs$S0,
        K = K,
        r = mtm_var_inputs$r,
        q = mtm_var_inputs$q,
        mid_vol = mtm_var_inputs$mid_vol,
        T = mtm_var_inputs$T,
        type = type,
        var_mu = mtm_var_inputs$var_mu,
        var_sigma = mtm_var_inputs$var_sigma,
        horizon_days = mtm_var_inputs$horizon_days,
        trading_days = mtm_var_inputs$trading_days,
        n_scenarios = mtm_var_inputs$n_scenarios,
        return_model = mtm_var_inputs$return_model,
        terminal_floor = mtm_var_inputs$terminal_floor,
        conf_levels = mtm_var_inputs$conf_levels
      )
    }
  }
  mtm_var_names <- grep("^(var|es)_", names(mtm_var), value = TRUE)
  mtm_var_values <- if (length(mtm_var_names) > 0) {
    unlist(mtm_var[mtm_var_names], use.names = FALSE)
  } else {
    numeric(0)
  }
  mtm_var_constraints <- normalize_mtm_loss_constraints(
    mtm_var_constraints,
    "mtm_var_constraints"
  )
  mtm_es_constraints <- normalize_mtm_loss_constraints(
    mtm_es_constraints,
    "mtm_es_constraints"
  )
  mtm_var_violation <- numeric(0)
  mtm_var_violation_names <- character(0)
  mtm_es_violation <- numeric(0)
  mtm_es_violation_names <- character(0)
  max_mtm_var_violation <- 0.0
  max_mtm_es_violation <- 0.0
  if (nrow(mtm_var_constraints) > 0) {
    mtm_var_constraint_names <- mtm_constraint_col_names("var", mtm_var_constraints$conf_level)
    mtm_var_violation_names <- mtm_constraint_col_names(
      "mtm_var_violation",
      mtm_var_constraints$conf_level
    )
    mtm_var_loss <- vapply(mtm_var_constraint_names, function(nm) mtm_var[[nm]], numeric(1))
    mtm_var_violation <- pmax(mtm_var_loss - mtm_var_constraints$max_loss, 0.0)
    max_mtm_var_violation <- max(mtm_var_violation)
  }
  if (nrow(mtm_es_constraints) > 0) {
    mtm_es_constraint_names <- mtm_constraint_col_names("es", mtm_es_constraints$conf_level)
    mtm_es_violation_names <- mtm_constraint_col_names(
      "mtm_es_violation",
      mtm_es_constraints$conf_level
    )
    mtm_es_loss <- vapply(mtm_es_constraint_names, function(nm) mtm_var[[nm]], numeric(1))
    mtm_es_violation <- pmax(mtm_es_loss - mtm_es_constraints$max_loss, 0.0)
    max_mtm_es_violation <- max(mtm_es_violation)
  }
  hedged_stats <- list(
    hedged_mean_wealth = NA_real_,
    hedged_sd_wealth = NA_real_,
    hedged_sharpe = NA_real_,
    hedged_skew = NA_real_,
    hedged_ex.kurt = NA_real_,
    hedged_mean_pnl = NA_real_,
    hedged_sd_pnl = NA_real_,
    hedged_avg_abs_final_delta = NA_real_
  )
  if (!is.null(delta_hedge_inputs) && delta_hedge_inputs$steps > 0L) {
    hedged_stats <- delta_hedged_portfolio_stats(
      contracts = contracts,
      cash = cash_after_trade,
      budget = budget,
      S0 = delta_hedge_inputs$S0,
      K = K,
      r = delta_hedge_inputs$r,
      q = delta_hedge_inputs$q,
      mid_vol = delta_hedge_inputs$mid_vol,
      T = delta_hedge_inputs$T,
      type = type,
      hedge_steps = delta_hedge_inputs$steps,
      hedge_paths = delta_hedge_inputs$paths,
      hedge_mu = delta_hedge_inputs$mu,
      hedge_sigma = delta_hedge_inputs$sigma,
      hedge_seed = delta_hedge_inputs$seed,
      stock_transaction_cost = delta_hedge_inputs$stock_transaction_cost
    )
  }

  list(
    invested_weight = sum(position_weight),
    trade_contracts = trade_contracts,
    trade_cost = trade_cost,
    cash_after_trade = cash_after_trade,
    cash_fraction = cash_after_trade / budget,
    mean_wealth = mean_wealth,
    variance_wealth = variance_wealth,
    sd_wealth = sd_wealth,
    edge_wealth = edge_wealth,
    sharpe = sharpe,
    adj_sharpe = as.numeric(adj_sharpe),
    skew = as.numeric(payoff_moments["skew"]),
    ex.kurt = as.numeric(payoff_moments["ex.kurt"]),
    min_wealth = min_wealth,
    min_wealth_violation = min_wealth_violation,
    tail_slope = tail_slope_value,
    max_es_violation = max_es_violation,
    mean_variance_objective = edge_wealth / budget -
      risk_aversion * variance_wealth / budget^2,
    gross_position_weight = gross_position_weight,
    dependency_objective = dependency_objective,
    dependency_penalty = dependency_penalty,
    portfolio_delta = portfolio_delta,
    portfolio_gamma = portfolio_gamma,
    portfolio_vega = portfolio_vega,
    delta_dollars = delta_dollars,
    gamma_dollars_1pct = gamma_dollars_1pct,
    vega_1pct = vega_1pct,
    delta_bound_violation = delta_bound_violation,
    current_mtm_wealth = mtm_var$current_mtm_wealth,
    mean_mtm_pnl = mtm_var$mean_mtm_pnl,
    sd_mtm_pnl = mtm_var$sd_mtm_pnl,
    mtm_var_names = mtm_var_names,
    mtm_var_values = mtm_var_values,
    max_mtm_var_violation = max_mtm_var_violation,
    max_mtm_es_violation = max_mtm_es_violation,
    mtm_var_violation = mtm_var_violation,
    mtm_var_violation_names = mtm_var_violation_names,
    mtm_es_violation = mtm_es_violation,
    mtm_es_violation_names = mtm_es_violation_names,
    hedged_mean_wealth = hedged_stats$hedged_mean_wealth,
    hedged_sd_wealth = hedged_stats$hedged_sd_wealth,
    hedged_sharpe = hedged_stats$hedged_sharpe,
    hedged_skew = hedged_stats$hedged_skew,
    hedged_ex.kurt = hedged_stats$hedged_ex.kurt,
    hedged_mean_pnl = hedged_stats$hedged_mean_pnl,
    hedged_sd_pnl = hedged_stats$hedged_sd_pnl,
    hedged_avg_abs_final_delta = hedged_stats$hedged_avg_abs_final_delta,
    port_corr = port_corr,
    trade_price = trade_price,
    payoff_per_dollar = payoff_per_dollar,
    edge = edge,
    instrument_sharpe = instrument_sharpe,
    instrument_adj_sharpe = instrument_adj_sharpe,
    trade_cost_by_instrument = trade_cost_by_instrument,
    trade_weight = trade_weight,
    position_value = position_value,
    position_weight = position_weight,
    es_loss = es_loss,
    es_names = es_names,
    es_violation = es_violation,
    es_violation_names = es_violation_names
  )
}

portfolio_objective_value <- function(optimization_objective,
                                      portfolio_eval,
                                      m,
                                      v,
                                      K,
                                      type,
                                      prices,
                                      budget,
                                      rf_growth,
                                      contracts,
                                      utility_gamma = NA_real_,
                                      initial_contracts = rep(0.0, length(contracts)),
                                      bid_prices = prices,
                                      ask_prices = prices,
                                      payoff_scenarios = NULL) {
  if (optimization_objective == "mean_variance") {
    return(list(
      name = "mean_variance_objective",
      value = portfolio_eval$mean_variance_objective
    ))
  }

  if (optimization_objective == "sharpe") {
    return(list(
      name = "sharpe_objective",
      value = portfolio_eval$sharpe
    ))
  }

  if (optimization_objective == "adjusted_sharpe" ||
      optimization_objective == "sharpe_adj") {
    return(list(
      name = "adjusted_sharpe_objective",
      value = portfolio_eval$adj_sharpe
    ))
  }

  if (optimization_objective == "expected_utility") {
    return(list(
      name = "expected_utility_objective",
      value = expected_utility_for_contracts(
        m = m,
        v = v,
        K = K,
        type = type,
        prices = prices,
        budget = budget,
        rf_growth = rf_growth,
        contracts = contracts,
        gamma = utility_gamma,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        payoff_scenarios = payoff_scenarios
      )
    ))
  }

  if (optimization_objective == "dependency_penalty") {
    return(list(
      name = "dependency_objective",
      value = portfolio_eval$dependency_objective
    ))
  }

  stop("optimization_objective must be 'mean_variance', 'sharpe', 'adjusted_sharpe', 'sharpe_adj', 'expected_utility', or 'dependency_penalty'")
}

option_instrument_index <- function(K, type, strike, option_type) {
  idx <- which(K == strike & type == option_type)
  if (length(idx) != 1L) {
    stop("Expected exactly one ", option_type, " at strike ", strike)
  }
  idx
}

option_structure_contract_vector <- function(K,
                                             type,
                                             structure_type,
                                             strikes,
                                             units = 1.0) {
  contracts <- numeric(length(K))
  add_leg <- function(option_type, strike, qty) {
    idx <- option_instrument_index(K, type, strike, option_type)
    contracts[idx] <<- contracts[idx] + units * qty
  }

  if (structure_type == "straddle") {
    add_leg("call", strikes[1], 1.0)
    add_leg("put", strikes[1], 1.0)
  } else if (structure_type == "strangle") {
    add_leg("put", strikes[1], 1.0)
    add_leg("call", strikes[2], 1.0)
  } else if (structure_type == "call_spread") {
    add_leg("call", strikes[1], 1.0)
    add_leg("call", strikes[2], -1.0)
  } else if (structure_type == "put_spread") {
    add_leg("put", strikes[1], -1.0)
    add_leg("put", strikes[2], 1.0)
  } else if (structure_type == "butterfly") {
    add_leg("call", strikes[1], 1.0)
    add_leg("call", strikes[2], -2.0)
    add_leg("call", strikes[3], 1.0)
  } else if (structure_type == "iron_condor") {
    add_leg("put", strikes[1], 1.0)
    add_leg("put", strikes[2], -1.0)
    add_leg("call", strikes[3], -1.0)
    add_leg("call", strikes[4], 1.0)
  } else {
    stop("Unsupported structure_type: ", structure_type)
  }

  contracts
}

generate_option_structure_candidates <- function(K,
                                                 type,
                                                 structure_types = c("straddle", "butterfly", "iron_condor"),
                                                 units = 1.0,
                                                 require_even_butterfly = TRUE) {
  strikes <- sort(unique(K[K > 0]))
  out <- list()
  add_candidate <- function(structure_type, strikes_for_structure) {
    contracts <- option_structure_contract_vector(
      K = K,
      type = type,
      structure_type = structure_type,
      strikes = strikes_for_structure,
      units = units
    )
    out[[length(out) + 1L]] <<- list(
      structure_type = structure_type,
      strikes = paste(strikes_for_structure, collapse = "/"),
      contracts = contracts
    )
  }

  for (structure_type in structure_types) {
    if (structure_type == "straddle") {
      for (k in strikes) add_candidate(structure_type, k)
    } else if (structure_type == "strangle") {
      if (length(strikes) >= 2L) {
        for (i in seq_len(length(strikes) - 1L)) {
          for (j in (i + 1L):length(strikes)) {
            add_candidate(structure_type, c(strikes[i], strikes[j]))
          }
        }
      }
    } else if (structure_type %in% c("call_spread", "put_spread")) {
      if (length(strikes) >= 2L) {
        for (i in seq_len(length(strikes) - 1L)) {
          for (j in (i + 1L):length(strikes)) {
            add_candidate(structure_type, c(strikes[i], strikes[j]))
          }
        }
      }
    } else if (structure_type == "butterfly") {
      if (length(strikes) >= 3L) {
        combos <- combn(strikes, 3)
        for (j in seq_len(ncol(combos))) {
          ks <- combos[, j]
          if (!require_even_butterfly || abs((ks[2] - ks[1]) - (ks[3] - ks[2])) < 1e-8) {
            add_candidate(structure_type, ks)
          }
        }
      }
    } else if (structure_type == "iron_condor") {
      if (length(strikes) >= 4L) {
        combos <- combn(strikes, 4)
        for (j in seq_len(ncol(combos))) {
          add_candidate(structure_type, combos[, j])
        }
      }
    } else {
      stop("Unsupported structure_type: ", structure_type)
    }
  }

  out
}

evaluate_option_structures <- function(candidates,
                                       m,
                                       v,
                                       K,
                                       type,
                                       option_price,
                                       bid_price,
                                       ask_price,
                                       mid_price,
                                       vol_quotes,
                                       greeks,
                                       S0,
                                       r,
                                       q,
                                       T,
                                       expected_payoff,
                                       cov_payoff,
                                       dependency_matrix = NULL,
                                       payoff_stats,
                                       payoff_scenarios,
                                       payoff_moment_scenarios = NULL,
                                       payoff_grid,
                                       tail_slope,
                                       cvar_constraints = NULL,
                                       cvar_payoff_scenarios = NULL,
                                       optimization_objective = "sharpe",
                                       utility_gamma = NA_real_,
                                       budget = 1000.0,
                                       rf_growth,
                                       risk_aversion = 0.0,
                                       min_terminal_wealth = -Inf,
                                       initial_contracts = rep(0.0, length(K)),
                                       portfolio_delta_bounds = c(NA_real_, NA_real_),
                                       mtm_var_inputs = NULL,
                                       mtm_var_constraints = NULL,
                                       mtm_es_constraints = NULL,
                                       delta_hedge_inputs = NULL,
                                       top_n = 10L,
                                       print_zero_weight_options = FALSE,
                                       zero_weight_tol = 1e-12) {
  rows <- data.frame()
  runs <- list()

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    contracts <- candidate$contracts
    portfolio_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = contracts,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = delta_hedge_inputs
    )
    objective_eval <- portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = portfolio_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = option_price,
      budget = budget,
      rf_growth = rf_growth,
      contracts = contracts,
      utility_gamma = utility_gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      payoff_scenarios = payoff_scenarios
    )
    table_result <- option_portfolio_table(
      K = K,
      type = type,
      vol_quotes = vol_quotes,
      mid_price = mid_price,
      bid_price = bid_price,
      ask_price = ask_price,
      greeks = greeks,
      initial_contracts = initial_contracts,
      contracts = contracts,
      portfolio_eval = portfolio_eval,
      print_zero_weight_options = print_zero_weight_options,
      zero_weight_tol = zero_weight_tol
    )

    row <- data.frame(
      structure_id = i,
      structure_type = candidate$structure_type,
      strikes = candidate$strikes,
      objective = optimization_objective,
      obj_value = objective_eval$value,
      trade_cost = portfolio_eval$trade_cost,
      cash_after_trade = portfolio_eval$cash_after_trade,
      mean_wealth = portfolio_eval$mean_wealth,
      sd_wealth = portfolio_eval$sd_wealth,
      sharpe = portfolio_eval$sharpe,
      adj_sharpe = portfolio_eval$adj_sharpe,
      skew = portfolio_eval$skew,
      ex.kurt = portfolio_eval$ex.kurt,
      min_wealth = portfolio_eval$min_wealth,
      max_es_violation = portfolio_eval$max_es_violation,
      portfolio_delta = portfolio_eval$portfolio_delta,
      portfolio_gamma = portfolio_eval$portfolio_gamma,
      portfolio_vega = portfolio_eval$portfolio_vega,
      current_mtm_wealth = portfolio_eval$current_mtm_wealth,
      mean_mtm_pnl = portfolio_eval$mean_mtm_pnl,
      sd_mtm_pnl = portfolio_eval$sd_mtm_pnl,
      hedged_mean_wealth = portfolio_eval$hedged_mean_wealth,
      hedged_sd_wealth = portfolio_eval$hedged_sd_wealth,
      hedged_sharpe = portfolio_eval$hedged_sharpe
    )
    for (j in seq_along(portfolio_eval$mtm_var_values)) {
      row[[portfolio_eval$mtm_var_names[j]]] <- portfolio_eval$mtm_var_values[j]
    }
    rows <- rbind(rows, row)
    runs[[length(runs) + 1L]] <- list(
      candidate = candidate,
      contracts = contracts,
      portfolio_eval = portfolio_eval,
      objective_eval = objective_eval,
      table_raw = table_result$raw,
      table_print = table_result$print
    )
  }

  rows <- rows[order(rows$obj_value, decreasing = TRUE), , drop = FALSE]
  rownames(rows) <- NULL
  top_n <- min(as.integer(top_n), nrow(rows))
  top_rows <- rows[seq_len(top_n), , drop = FALSE]
  list(
    summary_rows = rows,
    summary_print = format_option_optimization_summary(rows),
    top_summary = format_option_optimization_summary(top_rows),
    runs = runs,
    best_run = runs[[rows$structure_id[1]]]
  )
}

active_option_indices <- function(K, type) {
  which(!(type == "call" & K == 0))
}

prune_option_contracts <- function(contracts,
                                   portfolio_eval,
                                   K,
                                   type,
                                   max_active_option_positions = NA_integer_,
                                   prune_positions_by = "position_weight",
                                   zero_weight_tol = 1e-12) {
  if (is.na(max_active_option_positions)) {
    return(list(
      contracts = contracts,
      active_option_positions = sum(abs(contracts[active_option_indices(K, type)]) > zero_weight_tol),
      pruned_option_positions = 0L
    ))
  }

  max_active_option_positions <- as.integer(max_active_option_positions)
  if (max_active_option_positions < 0L) {
    stop("max_active_option_positions must be nonnegative or NA")
  }
  if (!prune_positions_by %in% c("position_weight", "trade_weight", "contracts")) {
    stop("prune_positions_by must be 'position_weight', 'trade_weight', or 'contracts'")
  }

  option_idx <- active_option_indices(K, type)
  active_idx <- option_idx[abs(contracts[option_idx]) > zero_weight_tol]
  active_count <- length(active_idx)
  if (active_count <= max_active_option_positions) {
    return(list(
      contracts = contracts,
      active_option_positions = active_count,
      pruned_option_positions = 0L
    ))
  }

  size <- switch(
    prune_positions_by,
    position_weight = abs(portfolio_eval$position_weight),
    trade_weight = abs(portfolio_eval$trade_weight),
    contracts = abs(contracts)
  )
  keep_idx <- active_idx[
    order(size[active_idx], decreasing = TRUE)[seq_len(max_active_option_positions)]
  ]
  pruned_contracts <- contracts
  prune_idx <- setdiff(active_idx, keep_idx)
  pruned_contracts[prune_idx] <- 0.0

  list(
    contracts = pruned_contracts,
    active_option_positions = length(keep_idx),
    pruned_option_positions = length(prune_idx)
  )
}

integer_contract_candidates <- function(continuous_contracts,
                                        initial_contracts,
                                        neighborhood = 1L,
                                        max_search_instruments = 6L) {
  trade <- continuous_contracts - initial_contracts
  rounded_trade <- round(trade)
  candidates <- list(initial_contracts + rounded_trade)
  fractional <- abs(trade - rounded_trade)
  search_idx <- order(fractional, decreasing = TRUE)
  search_idx <- search_idx[fractional[search_idx] > 1e-8]
  search_idx <- head(search_idx, max(0L, as.integer(max_search_instruments)))

  if (length(search_idx) == 0L || neighborhood < 1L) {
    return(candidates)
  }

  offsets <- seq.int(-as.integer(neighborhood), as.integer(neighborhood))
  grid <- expand.grid(rep(list(offsets), length(search_idx)))
  names(grid) <- as.character(search_idx)
  for (i in seq_len(nrow(grid))) {
    candidate_trade <- rounded_trade
    candidate_trade[search_idx] <- rounded_trade[search_idx] + as.numeric(grid[i, ])
    candidates[[length(candidates) + 1L]] <- initial_contracts + candidate_trade
  }

  candidate_matrix <- do.call(rbind, candidates)
  candidate_matrix <- unique(candidate_matrix)
  split(candidate_matrix, row(candidate_matrix))
}

portfolio_constraints_feasible <- function(portfolio_eval,
                                           contracts,
                                           constraint_mode,
                                           budget,
                                           max_invested_weight,
                                           force_full_investment,
                                           enforce_mtm_constraints = FALSE,
                                           tol = 1e-8) {
  if (!is.finite(portfolio_eval$trade_cost) ||
      portfolio_eval$trade_cost > budget + tol) {
    return(FALSE)
  }

  if (max(abs(portfolio_eval$position_weight), na.rm = TRUE) >
      max_invested_weight + tol) {
    return(FALSE)
  }

  if (portfolio_eval$delta_bound_violation > tol) {
    return(FALSE)
  }

  if (constraint_mode == "long_only") {
    if (any(contracts < -tol)) {
      return(FALSE)
    }
    if (portfolio_eval$invested_weight > max_invested_weight + tol) {
      return(FALSE)
    }
  }

  if (constraint_mode == "long_short") {
    if (sum(abs(portfolio_eval$position_weight), na.rm = TRUE) >
        max_invested_weight + tol) {
      return(FALSE)
    }
  }

  if (constraint_mode == "nonnegative_terminal") {
    if (portfolio_eval$min_wealth_violation > tol ||
        portfolio_eval$tail_slope > tol ||
        portfolio_eval$max_es_violation > tol) {
      return(FALSE)
    }
  }

  if (enforce_mtm_constraints &&
      (portfolio_eval$max_mtm_var_violation > tol ||
       portfolio_eval$max_mtm_es_violation > tol)) {
    return(FALSE)
  }

  TRUE
}

integer_repair_contracts <- function(continuous_contracts,
                                     initial_contracts,
                                     constraint_mode,
                                     optimization_objective,
                                     utility_gamma,
                                     m,
                                     v,
                                     K,
                                     type,
                                     option_price,
                                     bid_price,
                                     ask_price,
                                     budget,
                                     rf_growth,
                                     expected_payoff,
                                     cov_payoff,
                                     payoff_stats,
                                     payoff_scenarios,
                                     payoff_moment_scenarios,
                                     payoff_grid,
                                     tail_slope,
                                     cvar_constraints,
                                     cvar_payoff_scenarios,
                                     min_terminal_wealth,
                                     risk_aversion,
                                     dependency_matrix,
                                     greeks,
                                     S0,
                                     portfolio_delta_bounds,
                                     mtm_var_inputs,
                                     mtm_var_constraints,
                                     mtm_es_constraints,
                                     delta_hedge_inputs,
                                     max_invested_weight,
                                     force_full_investment,
                                     neighborhood = 1L,
                                     max_search_instruments = 6L,
                                     tol = 1e-8) {
  if (any(abs(initial_contracts - round(initial_contracts)) > tol)) {
    stop("integer_contracts requires integer initial_contracts")
  }

  candidates <- integer_contract_candidates(
    continuous_contracts = continuous_contracts,
    initial_contracts = initial_contracts,
    neighborhood = neighborhood,
    max_search_instruments = max_search_instruments
  )

  best <- NULL
  best_eval <- NULL
  best_objective <- NULL
  feasible_count <- 0L
  for (candidate in candidates) {
    candidate <- as.numeric(candidate)
    if (any(abs(candidate - initial_contracts -
            round(candidate - initial_contracts)) > tol)) {
      next
    }
    portfolio_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = candidate,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = NULL
    )
    if (!portfolio_constraints_feasible(
      portfolio_eval = portfolio_eval,
      contracts = candidate,
      constraint_mode = constraint_mode,
      budget = budget,
      max_invested_weight = max_invested_weight,
      force_full_investment = force_full_investment,
      enforce_mtm_constraints = TRUE,
      tol = tol
    )) {
      next
    }

    objective_eval <- portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = portfolio_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = option_price,
      budget = budget,
      rf_growth = rf_growth,
      contracts = candidate,
      utility_gamma = utility_gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      payoff_scenarios = payoff_scenarios
    )
    feasible_count <- feasible_count + 1L
    if (is.null(best_objective) || objective_eval$value > best_objective$value) {
      best <- candidate
      best_eval <- portfolio_eval
      best_objective <- objective_eval
    }
  }

  if (is.null(best)) {
    best <- initial_contracts
    best_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = best,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = NULL
    )
    best_objective <- portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = best_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = option_price,
      budget = budget,
      rf_growth = rf_growth,
      contracts = best,
      utility_gamma = utility_gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      payoff_scenarios = payoff_scenarios
    )
  }

  if (!is.null(delta_hedge_inputs) && delta_hedge_inputs$steps > 0L) {
    best_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = best,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = delta_hedge_inputs
    )
  }

  list(
    contracts = best,
    portfolio_eval = best_eval,
    objective_eval = best_objective,
    candidate_count = length(candidates),
    feasible_count = feasible_count,
    used_initial_fallback = feasible_count == 0L
  )
}

option_portfolio_table <- function(K,
                                   type,
                                   underlying = NULL,
                                   expiry = NULL,
                                   vol_quotes,
                                   mid_price,
                                   bid_price,
                                   ask_price,
                                   greeks,
                                   initial_contracts,
                                   contracts,
                                   portfolio_eval,
                                   print_zero_weight_options = TRUE,
                                   zero_weight_tol = 1e-12) {
  if (is.null(greeks)) {
    greeks <- data.frame(
      delta = rep(NA_real_, length(K)),
      gamma = rep(NA_real_, length(K)),
      vega = rep(NA_real_, length(K))
    )
  }
  out <- data.frame(
    type = type,
    strike = K,
    mid_vol = round(vol_quotes$mid_vol, 4),
    bid_vol = round(vol_quotes$bid_vol, 4),
    ask_vol = round(vol_quotes$ask_vol, 4),
    mid_price = round(mid_price, 3),
    bid_price = round(bid_price, 3),
    ask_price = round(ask_price, 3),
    trade_price = round(portfolio_eval$trade_price, 3),
    payoff_per_dollar = round(portfolio_eval$payoff_per_dollar, 3),
    edge = round(portfolio_eval$edge, 3),
    sharpe = round(portfolio_eval$instrument_sharpe, 3),
    adj_sharpe = round(portfolio_eval$instrument_adj_sharpe, 3),
    delta = round(greeks$delta, 4),
    gamma = round(greeks$gamma, 6),
    vega = round(greeks$vega, 4),
    port_corr = round(portfolio_eval$port_corr, 3),
    initial_contracts = initial_contracts,
    trade_contracts = portfolio_eval$trade_contracts,
    trade_cost = portfolio_eval$trade_cost_by_instrument,
    trade_weight = portfolio_eval$trade_weight,
    position_value = portfolio_eval$position_value,
    position_weight = portfolio_eval$position_weight,
    contracts = contracts
  )
  leading_cols <- c("type", "strike")
  extra <- data.frame(row_id = seq_along(K))
  if (!is.null(underlying)) {
    extra$underlying <- underlying
  }
  if (!is.null(expiry)) {
    extra$expiry <- expiry
  }
  extra$row_id <- NULL
  if (ncol(extra) > 0L) {
    out <- cbind(
      out[, leading_cols, drop = FALSE],
      extra,
      out[, setdiff(names(out), leading_cols), drop = FALSE]
    )
  }

  print_out <- out
  if (!print_zero_weight_options) {
    print_out <- out[abs(out$position_weight) > zero_weight_tol, , drop = FALSE]
  }

  numeric_4dp_columns <- c(
    "initial_contracts",
    "trade_contracts",
    "trade_cost",
    "trade_weight",
    "position_value",
    "position_weight",
    "contracts"
  )
  for (col in numeric_4dp_columns) {
    print_out[[col]] <- sprintf("%.4f", print_out[[col]])
  }

  list(raw = out, print = print_out)
}

option_objective_runs <- function(optimization_objectives,
                                  risk_aversion_utility) {
  objective_runs <- data.frame(objective = character(), utility_gamma = numeric())
  for (obj in optimization_objectives) {
    if (obj == "expected_utility") {
      objective_runs <- rbind(
        objective_runs,
        data.frame(objective = rep(obj, length(risk_aversion_utility)),
                   utility_gamma = risk_aversion_utility)
      )
    } else {
      objective_runs <- rbind(
        objective_runs,
        data.frame(objective = obj, utility_gamma = NA_real_)
      )
    }
  }
  objective_runs
}

print_option_optimization_run <- function(run) {
  cat("Option solution\n")
  cat("constraint_mode:", run$constraint_mode, "\n")
  cat("optimization_objective:", run$objective, "\n")
  if (run$objective == "expected_utility") {
    cat("utility_gamma:", run$utility_gamma, "\n")
  }
  cat("rf_growth:", round(run$rf_growth, 6), "\n")
  print(run$table_print, row.names = FALSE)
  cat("\ninvested_weight:", round(run$portfolio_eval$invested_weight, 6), "\n")
  cat("gross_position_weight:", round(run$portfolio_eval$gross_position_weight, 6), "\n")
  cat("trade_cost:", round(run$portfolio_eval$trade_cost, 6), "\n")
  cat("cash_after_trade:", round(run$portfolio_eval$cash_after_trade, 6), "\n")
  cat("cash_fraction:", round(run$portfolio_eval$cash_fraction, 6), "\n")
  cat("mean_wealth:", round(run$portfolio_eval$mean_wealth, 6), "\n")
  cat("sd_wealth:", round(run$portfolio_eval$sd_wealth, 6), "\n")
  cat("sharpe:", round(run$portfolio_eval$sharpe, 6), "\n")
  cat("adj_sharpe:", round(run$portfolio_eval$adj_sharpe, 6), "\n")
  cat("skew:", round(run$portfolio_eval$skew, 6), "\n")
  cat("ex.kurt:", round(run$portfolio_eval$ex.kurt, 6), "\n")
  cat("min_wealth:", round(run$portfolio_eval$min_wealth, 6), "\n")
  cat("min_wealth_violation:", round(run$portfolio_eval$min_wealth_violation, 6), "\n")
  cat("tail_slope:", round(run$portfolio_eval$tail_slope, 6), "\n")
  if (!is.na(run$portfolio_eval$portfolio_delta)) {
    cat("portfolio_delta:", round(run$portfolio_eval$portfolio_delta, 6), "\n")
    cat("portfolio_gamma:", round(run$portfolio_eval$portfolio_gamma, 6), "\n")
    cat("portfolio_vega:", round(run$portfolio_eval$portfolio_vega, 6), "\n")
    cat("delta_dollars:", round(run$portfolio_eval$delta_dollars, 6), "\n")
    cat("gamma_dollars_1pct:", round(run$portfolio_eval$gamma_dollars_1pct, 6), "\n")
    cat("vega_1pct:", round(run$portfolio_eval$vega_1pct, 6), "\n")
    cat("delta_bound_violation:", round(run$portfolio_eval$delta_bound_violation, 6), "\n")
  }
  cat("max_es_violation:", round(run$portfolio_eval$max_es_violation, 6), "\n")
  for (i in seq_along(run$portfolio_eval$es_loss)) {
    cat(run$portfolio_eval$es_names[i], ":", round(run$portfolio_eval$es_loss[i], 6), "\n")
  }
  for (i in seq_along(run$portfolio_eval$es_violation)) {
    cat(run$portfolio_eval$es_violation_names[i], ":", round(run$portfolio_eval$es_violation[i], 6), "\n")
  }
  if (!is.na(run$portfolio_eval$current_mtm_wealth)) {
    cat("current_mtm_wealth:", round(run$portfolio_eval$current_mtm_wealth, 6), "\n")
    cat("mean_mtm_pnl:", round(run$portfolio_eval$mean_mtm_pnl, 6), "\n")
    cat("sd_mtm_pnl:", round(run$portfolio_eval$sd_mtm_pnl, 6), "\n")
    for (i in seq_along(run$portfolio_eval$mtm_var_values)) {
      cat(run$portfolio_eval$mtm_var_names[i], ":",
          round(run$portfolio_eval$mtm_var_values[i], 6), "\n")
    }
    cat("max_mtm_var_violation:",
        round(run$portfolio_eval$max_mtm_var_violation, 6), "\n")
    cat("max_mtm_es_violation:",
        round(run$portfolio_eval$max_mtm_es_violation, 6), "\n")
    for (i in seq_along(run$portfolio_eval$mtm_var_violation)) {
      cat(run$portfolio_eval$mtm_var_violation_names[i], ":",
          round(run$portfolio_eval$mtm_var_violation[i], 6), "\n")
    }
    for (i in seq_along(run$portfolio_eval$mtm_es_violation)) {
      cat(run$portfolio_eval$mtm_es_violation_names[i], ":",
          round(run$portfolio_eval$mtm_es_violation[i], 6), "\n")
    }
  }
  if (!is.na(run$portfolio_eval$hedged_mean_wealth)) {
    cat("hedged_mean_wealth:", round(run$portfolio_eval$hedged_mean_wealth, 6), "\n")
    cat("hedged_sd_wealth:", round(run$portfolio_eval$hedged_sd_wealth, 6), "\n")
    cat("hedged_sharpe:", round(run$portfolio_eval$hedged_sharpe, 6), "\n")
    cat("hedged_skew:", round(run$portfolio_eval$hedged_skew, 6), "\n")
    cat("hedged_ex.kurt:", round(run$portfolio_eval$hedged_ex.kurt, 6), "\n")
    cat("hedged_mean_pnl:", round(run$portfolio_eval$hedged_mean_pnl, 6), "\n")
    cat("hedged_sd_pnl:", round(run$portfolio_eval$hedged_sd_pnl, 6), "\n")
    cat("hedged_avg_abs_final_delta:",
        round(run$portfolio_eval$hedged_avg_abs_final_delta, 6), "\n")
  }
  cat(run$obj_name, ":", round(run$obj_value, 6), "\n")
  if (run$obj_name != "mean_variance_objective") {
    cat("mean_variance_objective:",
        round(run$portfolio_eval$mean_variance_objective, 6), "\n")
  }
  if (!is.na(run$portfolio_eval$dependency_objective) &&
      run$obj_name != "dependency_objective") {
    cat("dependency_objective:",
        round(run$portfolio_eval$dependency_objective, 6), "\n")
    cat("dependency_penalty:",
        round(run$portfolio_eval$dependency_penalty, 6), "\n")
  } else if (!is.na(run$portfolio_eval$dependency_penalty)) {
    cat("dependency_penalty:",
        round(run$portfolio_eval$dependency_penalty, 6), "\n")
  }
  cat("continuous_obj_value:", round(run$continuous_obj_value, 6), "\n")
  cat("pruned_obj_value:", round(run$pruned_obj_value, 6), "\n")
  cat("prune_obj_delta:", round(run$prune_obj_delta, 6), "\n")
  cat("prune_obj_loss:", round(run$prune_obj_loss, 6), "\n")
  cat("active_option_positions:", run$active_option_positions, "\n")
  cat("pruned_option_positions:", run$pruned_option_positions, "\n")
  if (!is.na(run$integer_candidate_count)) {
    cat("integer_obj_delta:", round(run$integer_obj_delta, 6), "\n")
    cat("integer_obj_loss:", round(run$integer_obj_loss, 6), "\n")
    cat("integer_candidate_count:", run$integer_candidate_count, "\n")
    cat("integer_feasible_count:", run$integer_feasible_count, "\n")
    cat("integer_used_initial_fallback:", run$integer_used_initial_fallback, "\n")
  }
  cat("method_elapsed_seconds:", round(run$method_elapsed_seconds, 6), "\n")
  cat("\n")
}

select_option_optimization_run <- function(runs,
                                           constraint_mode = NA_character_,
                                           objective = NA_character_,
                                           utility_gamma = NA_real_,
                                           tol = 1e-10) {
  if (length(runs) == 0) {
    stop("No optimization runs are available")
  }

  for (run in runs) {
    constraint_ok <- is.na(constraint_mode) || run$constraint_mode == constraint_mode
    objective_ok <- is.na(objective) || run$objective == objective
    gamma_ok <- is.na(utility_gamma) ||
      (!is.na(run$utility_gamma) && abs(run$utility_gamma - utility_gamma) <= tol)
    if (constraint_ok && objective_ok && gamma_ok) {
      return(run)
    }
  }

  stop("No optimization run matched the requested plot selector")
}

plot_option_optimization_report <- function(file,
                                            selected_run,
                                            S0,
                                            r,
                                            q,
                                            T,
                                            K,
                                            type,
                                            market_vol,
                                            expected_payoff,
                                            rf_growth,
                                            width = 1200,
                                            height = 1000) {
  option_idx <- K > 0
  if (!any(option_idx)) {
    stop("No option rows are available to plot")
  }

  plot_K <- K[option_idx]
  plot_type <- type[option_idx]
  true_price <- exp(-r * T) * expected_payoff[option_idx]
  true_vol <- bs_option_implied_vol_vec(true_price, S0, plot_K, r, q, T, plot_type)
  market_plot_vol <- market_vol[option_idx]
  call_idx <- plot_type == "call"
  put_idx <- plot_type == "put"

  edge <- selected_run$portfolio_eval$edge[option_idx]
  table <- selected_run$table_raw
  bar_labels <- ifelse(table$type == "call" & table$strike == 0, "stock", paste(table$type, table$strike))
  bar_colors <- ifelse(table$type == "call", "steelblue3", "firebrick3")

  png(file, width = width, height = height)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1.1, 0.9, 1.1))
  par(mar = c(4.0, 4.5, 3.2, 1.2))

  ylim_vol <- range(c(market_plot_vol, true_vol), finite = TRUE)
  plot(
    plot_K,
    market_plot_vol,
    type = "n",
    ylim = ylim_vol,
    xlab = "Strike",
    ylab = "Implied vol",
    main = "Market vs model-implied true volatility"
  )
  grid(col = "gray88")
  lines(plot_K[call_idx], market_plot_vol[call_idx], type = "b", pch = 16, col = "steelblue3", lty = 2)
  lines(plot_K[put_idx], market_plot_vol[put_idx], type = "b", pch = 17, col = "firebrick3", lty = 2)
  lines(plot_K[call_idx], true_vol[call_idx], type = "b", pch = 1, col = "steelblue4", lwd = 2)
  lines(plot_K[put_idx], true_vol[put_idx], type = "b", pch = 2, col = "firebrick4", lwd = 2)
  legend(
    "topright",
    legend = c("market call", "market put", "true call", "true put"),
    col = c("steelblue3", "firebrick3", "steelblue4", "firebrick4"),
    lty = c(2, 2, 1, 1),
    pch = c(16, 17, 1, 2),
    bty = "n",
    cex = 0.9
  )

  ylim_edge <- range(edge, finite = TRUE)
  plot(
    plot_K,
    edge,
    type = "n",
    ylim = ylim_edge,
    xlab = "Strike",
    ylab = "Payoff edge",
    main = "True payoff edge versus cash growth"
  )
  abline(h = 0, col = "gray60", lty = 2)
  grid(col = "gray88")
  lines(plot_K[call_idx], edge[call_idx], type = "b", pch = 16, col = "steelblue3")
  lines(plot_K[put_idx], edge[put_idx], type = "b", pch = 17, col = "firebrick3")
  legend("topright", legend = c("call", "put"), col = c("steelblue3", "firebrick3"), pch = c(16, 17), lty = 1, bty = "n")

  par(mar = c(6.5, 4.5, 3.2, 1.2))
  bp <- barplot(
    table$position_weight,
    names.arg = bar_labels,
    col = bar_colors,
    border = NA,
    las = 2,
    ylab = "Position weight",
    main = paste(
      "Recommended positions:",
      selected_run$constraint_mode,
      selected_run$objective,
      if (!is.na(selected_run$utility_gamma)) paste0("gamma=", selected_run$utility_gamma) else ""
    )
  )
  abline(h = 0, col = "gray35")
  grid(nx = NA, ny = NULL, col = "gray88")
  legend("topright", legend = c("call/stock", "put"), fill = c("steelblue3", "firebrick3"), bty = "n")
  stats_text <- sprintf(
    "mean %.1f | sd %.1f | Sharpe %.3f | adj %.3f | min %.1f | ES1%% %.1f | cash %.3f",
    selected_run$portfolio_eval$mean_wealth,
    selected_run$portfolio_eval$sd_wealth,
    selected_run$portfolio_eval$sharpe,
    selected_run$portfolio_eval$adj_sharpe,
    selected_run$portfolio_eval$min_wealth,
    if (length(selected_run$portfolio_eval$es_loss) > 0) selected_run$portfolio_eval$es_loss[1] else NA_real_,
    selected_run$portfolio_eval$cash_fraction
  )
  mtext(stats_text, side = 1, line = 5.2, cex = 0.85)

  invisible(file)
}

format_option_optimization_summary <- function(summary_rows) {
  summary_print <- summary_rows
  num_cols <- vapply(summary_print, is.numeric, logical(1))
  summary_print[num_cols] <- lapply(summary_print[num_cols], round, 6)
  wealth_cols <- c(
    "mean_wealth",
    "sd_wealth",
    "sharpe",
    "adj_sharpe",
    "skew",
    "ex.kurt",
    "min_wealth",
    "min_wealth_violation",
    "tail_slope",
    "max_es_violation",
    "current_mtm_wealth",
    "mean_mtm_pnl",
    "sd_mtm_pnl",
    "max_mtm_var_violation",
    "max_mtm_es_violation",
    "hedged_mean_wealth",
    "hedged_sd_wealth",
    "hedged_sharpe",
    "hedged_skew",
    "hedged_ex.kurt",
    "hedged_mean_pnl",
    "hedged_sd_pnl",
    "hedged_avg_abs_final_delta",
    grep("^var_", names(summary_rows), value = TRUE),
    grep("^es_", names(summary_rows), value = TRUE),
    grep("^mtm_var_violation_", names(summary_rows), value = TRUE),
    grep("^mtm_es_violation_", names(summary_rows), value = TRUE),
    grep("^es_violation_", names(summary_rows), value = TRUE),
    grep("^es_loss_", names(summary_rows), value = TRUE)
  )
  for (col in wealth_cols) {
    if (col %in% names(summary_print)) {
      summary_print[[col]] <- round(summary_rows[[col]], 3)
    }
  }
  summary_print
}

option_optimizer_diagnostics <- function(opt) {
  diag <- opt$optimizer_diagnostics
  if (is.null(diag)) {
    diag <- list(
      optimizer = NA_character_,
      optimizer_status = NA_character_,
      optimizer_status_code = NA_integer_,
      optimizer_evals = NA_integer_,
      optimizer_fallback_used = NA
    )
  }
  diag
}

decode_nloptr_status <- function(status) {
  status <- as.integer(status)
  labels <- c(
    "1" = "success",
    "2" = "stopval_reached",
    "3" = "ftol_reached",
    "4" = "xtol_reached",
    "5" = "maxeval_reached",
    "6" = "maxtime_reached",
    "-1" = "failure",
    "-2" = "invalid_args",
    "-3" = "out_of_memory",
    "-4" = "roundoff_limited",
    "-5" = "forced_stop"
  )
  out <- labels[as.character(status)]
  if (is.na(out)) paste0("status_", status) else unname(out)
}

decode_optim_status <- function(status) {
  status <- as.integer(status)
  if (status == 0L) {
    "converged"
  } else if (status == 1L) {
    "maxit_reached"
  } else {
    paste0("convergence_", status)
  }
}

add_contract_start <- function(starts, x0, ui, ci, tol = 1e-8) {
  if (is.null(x0) || length(x0) != ncol(ui) || any(!is.finite(x0))) {
    return(starts)
  }
  if (!all(as.numeric(ui %*% x0 - ci) >= -tol)) {
    return(starts)
  }
  is_duplicate <- any(vapply(starts, function(existing) {
    length(existing) == length(x0) && max(abs(existing - x0)) <= tol
  }, logical(1)))
  if (!is_duplicate) {
    starts[[length(starts) + 1]] <- x0
  }
  starts
}

run_option_optimization_grid <- function(m,
                                         v,
                                         K,
                                         type,
                                         option_price,
                                         bid_price,
                                         ask_price,
                                         mid_price,
                                         vol_quotes,
                                         greeks,
                                         S0 = NA_real_,
                                         r = NA_real_,
                                         q = NA_real_,
                                         T = NA_real_,
                                         instrument_underlying = NULL,
                                         instrument_expiry = NULL,
                                         portfolio_delta_bounds = c(NA_real_, NA_real_),
                                         expected_payoff,
                                         cov_payoff,
                                         dependency_matrix = NULL,
                                         dependency_diagnostics = NULL,
                                         payoff_stats,
                                         payoff_scenarios,
                                         payoff_moment_scenarios = NULL,
                                         payoff_grid,
                                         tail_slope,
                                         lower_tail_slope = NULL,
                                         cvar_constraints = NULL,
                                         cvar_scenarios = NULL,
                                         cvar_payoff_scenarios = NULL,
                                         optimization_objectives,
                                         risk_aversion_utility,
                                         constraint_modes,
                                         budget,
                                         rf_growth,
                                         risk_aversion,
                                         min_terminal_wealth,
                                         initial_contracts,
                                         max_invested_weight = 1.0,
                                         force_full_investment = TRUE,
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
                                         mtm_var_inputs_override = NULL,
                                         mtm_var_constraints = NULL,
                                         mtm_es_constraints = NULL,
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
                                         tol = 1e-10,
                                         print_runs = TRUE) {
  payoff_per_dollar <- expected_payoff / option_price
  edge <- payoff_per_dollar - rf_growth
  cov_per_dollar <- cov_payoff / outer(option_price, option_price)
  portfolio_delta_bounds <- normalize_portfolio_delta_bounds(portfolio_delta_bounds)
  delta_for_constraints <- if (!is.null(greeks)) greeks$delta else NULL
  if (is.na(var_return_mu)) {
    var_return_mu <- if (is.finite(r) && is.finite(q)) r - q else 0.0
  }
  if (is.na(var_return_sigma)) {
    var_return_sigma <- if (is.finite(T) && T > 0.0) sqrt(v / T) else sqrt(v)
  }
  mtm_var_constraints <- normalize_mtm_loss_constraints(
    mtm_var_constraints,
    "mtm_var_constraints"
  )
  mtm_es_constraints <- normalize_mtm_loss_constraints(
    mtm_es_constraints,
    "mtm_es_constraints"
  )
  if ((nrow(mtm_var_constraints) > 0 || nrow(mtm_es_constraints) > 0) &&
      !isTRUE(report_mtm_var)) {
    stop("MTM VaR/ES constraints require report_mtm_var = TRUE")
  }
  var_conf_levels <- sort(unique(c(
    var_conf_levels,
    mtm_var_constraints$conf_level,
    mtm_es_constraints$conf_level
  )))
  mtm_var_inputs <- mtm_var_inputs_override
  if (report_mtm_var && is.null(mtm_var_inputs)) {
    mtm_var_inputs <- list(
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
  }
  delta_hedge_steps <- as.integer(delta_hedge_steps)
  delta_hedge_paths <- as.integer(delta_hedge_paths)
  if (delta_hedge_steps < 0L) {
    stop("delta_hedge_steps must be nonnegative")
  }
  if (is.na(delta_hedge_mu)) {
    delta_hedge_mu <- var_return_mu
  }
  if (is.na(delta_hedge_sigma)) {
    delta_hedge_sigma <- var_return_sigma
  }
  delta_hedge_inputs <- NULL
  if (delta_hedge_steps > 0L) {
    delta_hedge_inputs <- list(
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
  }
  if (is.null(dependency_matrix) &&
      "dependency_penalty" %in% optimization_objectives) {
    dependency_repair <- psd_repair_matrix(
      option_payoff_dependency_matrix(payoff_scenarios),
      eigen_floor = 1e-8
    )
    dependency_matrix <- dependency_repair$matrix
    dependency_diagnostics <- dependency_repair
  }
  dependency_diag <- if (is.null(dependency_diagnostics)) {
    list(
      raw_min_eigen = NA_real_,
      raw_max_eigen = NA_real_,
      repaired_min_eigen = NA_real_,
      repaired_max_eigen = NA_real_,
      condition_number = NA_real_,
      repaired = NA
    )
  } else {
    dependency_diagnostics
  }

  summary_rows <- data.frame()
  runs <- list()
  objective_runs <- option_objective_runs(
    optimization_objectives = optimization_objectives,
    risk_aversion_utility = risk_aversion_utility
  )

  for (constraint_mode in constraint_modes) {
  nonnegative_warm_starts <- list()
  for (run_idx in seq_len(nrow(objective_runs))) {
    optimization_objective <- objective_runs$objective[run_idx]
    utility_gamma <- objective_runs$utility_gamma[run_idx]
    method_start_time <- proc.time()

    if (constraint_mode == "nonnegative_terminal") {
      if (optimization_objective == "expected_utility") {
        opt <- expected_utility_option_opt(
          m = m,
          v = v,
          K = K,
          type = type,
          prices = option_price,
          bid_prices = bid_price,
          ask_prices = ask_price,
          budget = budget,
          rf_growth = rf_growth,
          gamma = utility_gamma,
          min_terminal_wealth = min_terminal_wealth,
          initial_contracts = initial_contracts,
          cvar_constraints = cvar_constraints,
          cvar_scenarios = cvar_scenarios,
          payoff_scenarios = payoff_scenarios,
          payoff_grid = payoff_grid,
          lower_tail_slope = lower_tail_slope,
          delta = delta_for_constraints,
          portfolio_delta_bounds = portfolio_delta_bounds,
          max_abs_position_weight = max_invested_weight,
          max_starts = nonnegative_expected_utility_max_starts,
          max_iter = nonnegative_expected_utility_max_iter,
          constrained_optimizer = constrained_optimizer,
          extra_starts = if (constrained_optimizer_use_warm_starts) nonnegative_warm_starts else list(),
          tol = tol
        )
        contracts <- opt$contracts
      } else {
        opt <- nonnegative_terminal_option_opt(
          m = m,
          v = v,
          K = K,
          type = type,
          prices = option_price,
          bid_prices = bid_price,
          ask_prices = ask_price,
          expected_payoff = expected_payoff,
          cov_payoff = cov_payoff,
          dependency_matrix = dependency_matrix,
          budget = budget,
          rf_growth = rf_growth,
          objective = optimization_objective,
          risk_aversion = risk_aversion,
          min_terminal_wealth = min_terminal_wealth,
          initial_contracts = initial_contracts,
          cvar_constraints = cvar_constraints,
          cvar_scenarios = cvar_scenarios,
          payoff_scenarios = payoff_moment_scenarios,
          payoff_grid = payoff_grid,
          lower_tail_slope = lower_tail_slope,
          delta = delta_for_constraints,
          portfolio_delta_bounds = portfolio_delta_bounds,
          max_abs_position_weight = max_invested_weight,
          constrained_optimizer = constrained_optimizer,
          max_starts = constrained_optimizer_max_starts,
          max_iter = constrained_optimizer_max_iter,
          extra_starts = if (constrained_optimizer_use_warm_starts) nonnegative_warm_starts else list(),
          tol = tol
        )
        contracts <- opt$contracts
      }
    } else if (constraint_mode == "long_only" &&
               optimization_objective == "mean_variance") {
      opt <- long_only_mean_variance_opt(
        edge = edge,
        covmat = cov_per_dollar,
        risk_aversion = risk_aversion,
        max_weight = max_invested_weight,
        force_full_investment = force_full_investment,
        tol = tol
      )
      contracts <- budget * opt$weights / option_price
    } else if (constraint_mode == "long_only" &&
               optimization_objective == "dependency_penalty") {
      if (is.null(dependency_matrix)) {
        stop("dependency_penalty objective requires dependency_matrix")
      }
      opt <- long_only_mean_variance_opt(
        edge = edge,
        covmat = dependency_matrix,
        risk_aversion = risk_aversion,
        max_weight = max_invested_weight,
        force_full_investment = force_full_investment,
        tol = tol
      )
      contracts <- budget * opt$weights / option_price
    } else if (constraint_mode == "long_only" &&
               optimization_objective == "sharpe") {
      opt <- long_only_max_sharpe_opt(
        edge = edge,
        covmat = cov_per_dollar,
        max_weight = max_invested_weight,
        tol = tol
      )
      contracts <- budget * opt$weights / option_price
    } else if (constraint_mode == "long_only" &&
               (optimization_objective == "adjusted_sharpe" ||
                optimization_objective == "sharpe_adj")) {
      opt <- long_only_max_adjusted_sharpe_opt(
        m = m,
        v = v,
        K = K,
        type = type,
        prices = option_price,
        edge = edge,
        covmat = cov_per_dollar,
        payoff_scenarios = payoff_moment_scenarios,
        max_weight = max_invested_weight,
        tol = tol
      )
      contracts <- budget * opt$weights / option_price
    } else if (constraint_mode == "long_only" &&
               optimization_objective == "expected_utility") {
      opt <- long_only_expected_utility_opt(
        m = m,
        v = v,
        K = K,
        type = type,
        prices = option_price,
        budget = budget,
        rf_growth = rf_growth,
        gamma = utility_gamma,
        initial_contracts = initial_contracts,
        bid_prices = bid_price,
        ask_prices = ask_price,
        payoff_scenarios = payoff_scenarios,
        max_weight = max_invested_weight,
        force_full_investment = force_full_investment
      )
      contracts <- opt$contracts
    } else if (constraint_mode == "long_short") {
      if (optimization_objective == "dependency_penalty" &&
          is.null(dependency_matrix)) {
        stop("dependency_penalty objective requires dependency_matrix")
      }
      opt <- long_short_portfolio_opt(
        optimization_objective = optimization_objective,
        m = m,
        v = v,
        K = K,
        type = type,
        bid_prices = bid_price,
        ask_prices = ask_price,
        budget = budget,
        rf_growth = rf_growth,
        expected_payoff = expected_payoff,
        cov_payoff = cov_payoff,
        payoff_stats = payoff_stats,
        payoff_scenarios = payoff_scenarios,
        payoff_moment_scenarios = payoff_moment_scenarios,
        payoff_grid = payoff_grid,
        tail_slope = tail_slope,
        risk_aversion = risk_aversion,
        dependency_matrix = dependency_matrix,
        gamma = utility_gamma,
        initial_contracts = initial_contracts,
        max_gross_weight = max_invested_weight,
        force_full_investment = force_full_investment,
        tol = tol
      )
      contracts <- opt$contracts
    } else {
      stop("Unsupported constraint_mode/objective combination")
    }

    if (constraint_mode == "nonnegative_terminal") {
      contracts <- scale_contracts_to_terminal_constraints(
        contracts = contracts,
        prices = option_price,
        budget = budget,
        rf_growth = rf_growth,
        payoff_grid = payoff_grid,
        lower_tail_slope = lower_tail_slope,
        min_terminal_wealth = min_terminal_wealth,
        cvar_constraints = cvar_constraints,
        cvar_payoff_scenarios = cvar_payoff_scenarios,
        delta = delta_for_constraints,
        portfolio_delta_bounds = portfolio_delta_bounds,
        initial_contracts = initial_contracts,
        bid_prices = bid_price,
        ask_prices = ask_price,
        tol = 1e-8
      )
      if (constrained_optimizer_use_warm_starts) {
        warm_delta_constraints <- portfolio_delta_constraint_rows(
          delta_for_constraints,
          portfolio_delta_bounds
        )
        nonnegative_warm_starts <- add_contract_start(
          nonnegative_warm_starts,
          contracts,
          ui = rbind(
            payoff_grid - rf_growth * matrix(option_price, nrow(payoff_grid), length(option_price), byrow = TRUE),
            tail_slope,
            if (!is.null(lower_tail_slope)) -lower_tail_slope else NULL,
            warm_delta_constraints$ui,
            -option_price,
            option_price,
            -diag(option_price),
            diag(option_price)
          ),
          ci = c(
            rep(min_terminal_wealth - budget * rf_growth, nrow(payoff_grid)),
            0.0,
            if (!is.null(lower_tail_slope)) 0.0 else NULL,
            warm_delta_constraints$ci,
            -budget,
            -budget,
            rep(-max_invested_weight * budget, length(option_price)),
            rep(-max_invested_weight * budget, length(option_price))
          )
        )
      }
    } else if (!is.null(delta_for_constraints) &&
               portfolio_delta_bound_violation(sum(delta_for_constraints * contracts),
                                               portfolio_delta_bounds) > 1e-8) {
      contracts <- scale_contracts_to_terminal_constraints(
        contracts = contracts,
        prices = option_price,
        budget = budget,
        rf_growth = rf_growth,
        payoff_grid = payoff_grid,
        lower_tail_slope = NULL,
        min_terminal_wealth = -Inf,
        cvar_constraints = NULL,
        cvar_payoff_scenarios = NULL,
        delta = delta_for_constraints,
        portfolio_delta_bounds = portfolio_delta_bounds,
        initial_contracts = initial_contracts,
        bid_prices = bid_price,
        ask_prices = ask_price,
        tol = 1e-8
      )
    }

    continuous_contracts <- contracts
    continuous_portfolio_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = contracts,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = delta_hedge_inputs
    )
    continuous_objective_eval <- portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = continuous_portfolio_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = option_price,
      budget = budget,
      rf_growth = rf_growth,
      contracts = contracts,
      utility_gamma = utility_gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      payoff_scenarios = payoff_scenarios
    )
    prune_result <- prune_option_contracts(
      contracts = continuous_contracts,
      portfolio_eval = continuous_portfolio_eval,
      K = K,
      type = type,
      max_active_option_positions = max_active_option_positions,
      prune_positions_by = prune_positions_by,
      zero_weight_tol = zero_weight_tol
    )
    pruned_contracts <- prune_result$contracts
    if (prune_repair_constraints && prune_result$pruned_option_positions > 0L) {
      pruned_contracts <- scale_contracts_to_terminal_constraints(
        contracts = pruned_contracts,
        prices = option_price,
        budget = budget,
        rf_growth = rf_growth,
        payoff_grid = payoff_grid,
        lower_tail_slope = if (constraint_mode == "nonnegative_terminal") lower_tail_slope else NULL,
        min_terminal_wealth = if (constraint_mode == "nonnegative_terminal") min_terminal_wealth else -Inf,
        cvar_constraints = if (constraint_mode == "nonnegative_terminal") cvar_constraints else NULL,
        cvar_payoff_scenarios = if (constraint_mode == "nonnegative_terminal") cvar_payoff_scenarios else NULL,
        delta = delta_for_constraints,
        portfolio_delta_bounds = portfolio_delta_bounds,
        initial_contracts = initial_contracts,
        bid_prices = bid_price,
        ask_prices = ask_price,
        tol = 1e-8
      )
    }
    pruned_portfolio_eval <- evaluate_option_portfolio(
      m = m,
      v = v,
      K = K,
      type = type,
      contracts = pruned_contracts,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      budget = budget,
      rf_growth = rf_growth,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      payoff_stats = payoff_stats,
      payoff_scenarios = payoff_moment_scenarios,
      payoff_grid = payoff_grid,
      tail_slope = tail_slope,
      cvar_constraints = cvar_constraints,
      cvar_payoff_scenarios = cvar_payoff_scenarios,
      min_terminal_wealth = min_terminal_wealth,
      risk_aversion = risk_aversion,
      dependency_matrix = dependency_matrix,
      greeks = greeks,
      S0 = S0,
      portfolio_delta_bounds = portfolio_delta_bounds,
      mtm_var_inputs = mtm_var_inputs,
      mtm_var_constraints = mtm_var_constraints,
      mtm_es_constraints = mtm_es_constraints,
      delta_hedge_inputs = delta_hedge_inputs
    )
    pruned_objective_eval <- portfolio_objective_value(
      optimization_objective = optimization_objective,
      portfolio_eval = pruned_portfolio_eval,
      m = m,
      v = v,
      K = K,
      type = type,
      prices = option_price,
      budget = budget,
      rf_growth = rf_growth,
      contracts = pruned_contracts,
      utility_gamma = utility_gamma,
      initial_contracts = initial_contracts,
      bid_prices = bid_price,
      ask_prices = ask_price,
      payoff_scenarios = payoff_scenarios
    )
    active_option_positions <- prune_result$active_option_positions
    pruned_option_positions <- prune_result$pruned_option_positions
    prune_obj_delta <- pruned_objective_eval$value - continuous_objective_eval$value
    prune_obj_loss <- max(-prune_obj_delta, 0.0)

    contracts <- pruned_contracts
    portfolio_eval <- pruned_portfolio_eval
    objective_eval <- pruned_objective_eval
    integer_candidate_count <- NA_integer_
    integer_feasible_count <- NA_integer_
    integer_used_initial_fallback <- NA
    if (integer_contracts) {
      integer_result <- integer_repair_contracts(
        continuous_contracts = pruned_contracts,
        initial_contracts = initial_contracts,
        constraint_mode = constraint_mode,
        optimization_objective = optimization_objective,
        utility_gamma = utility_gamma,
        m = m,
        v = v,
        K = K,
        type = type,
        option_price = option_price,
        bid_price = bid_price,
        ask_price = ask_price,
        budget = budget,
        rf_growth = rf_growth,
        expected_payoff = expected_payoff,
        cov_payoff = cov_payoff,
        payoff_stats = payoff_stats,
        payoff_scenarios = payoff_scenarios,
        payoff_moment_scenarios = payoff_moment_scenarios,
        payoff_grid = payoff_grid,
        tail_slope = tail_slope,
        cvar_constraints = cvar_constraints,
        cvar_payoff_scenarios = cvar_payoff_scenarios,
        min_terminal_wealth = min_terminal_wealth,
        risk_aversion = risk_aversion,
        dependency_matrix = dependency_matrix,
        greeks = greeks,
        S0 = S0,
        portfolio_delta_bounds = portfolio_delta_bounds,
        mtm_var_inputs = mtm_var_inputs,
        mtm_var_constraints = mtm_var_constraints,
        mtm_es_constraints = mtm_es_constraints,
        delta_hedge_inputs = delta_hedge_inputs,
        max_invested_weight = max_invested_weight,
        force_full_investment = force_full_investment,
        neighborhood = integer_rounding_neighborhood,
        max_search_instruments = integer_max_search_instruments,
        tol = 1e-8
      )
      contracts <- integer_result$contracts
      portfolio_eval <- integer_result$portfolio_eval
      objective_eval <- integer_result$objective_eval
      integer_candidate_count <- integer_result$candidate_count
      integer_feasible_count <- integer_result$feasible_count
      integer_used_initial_fallback <- integer_result$used_initial_fallback
    } else {
      portfolio_eval <- pruned_portfolio_eval
      objective_eval <- pruned_objective_eval
    }
    integer_obj_delta <- objective_eval$value - continuous_objective_eval$value
    integer_obj_loss <- max(-integer_obj_delta, 0.0)
    table_result <- option_portfolio_table(
      K = K,
      type = type,
      underlying = instrument_underlying,
      expiry = instrument_expiry,
      vol_quotes = vol_quotes,
      mid_price = mid_price,
      bid_price = bid_price,
      ask_price = ask_price,
      greeks = greeks,
      initial_contracts = initial_contracts,
      contracts = contracts,
      portfolio_eval = portfolio_eval,
      print_zero_weight_options = print_zero_weight_options,
      zero_weight_tol = zero_weight_tol
    )

    method_elapsed <- proc.time() - method_start_time
    method_elapsed_seconds <- as.numeric(method_elapsed["elapsed"])
    opt_diag <- option_optimizer_diagnostics(opt)
    if (is.na(opt_diag$optimizer)) {
      opt_diag$optimizer <- if (constraint_mode == "long_only") "long_only_internal" else "internal"
      opt_diag$optimizer_status <- "ok"
      opt_diag$optimizer_status_code <- NA_integer_
      opt_diag$optimizer_evals <- NA_integer_
      opt_diag$optimizer_fallback_used <- FALSE
    }
    run <- list(
      constraint_mode = constraint_mode,
      objective = optimization_objective,
      utility_gamma = utility_gamma,
      rf_growth = rf_growth,
      contracts = contracts,
      continuous_contracts = continuous_contracts,
      optimizer_result = opt,
      portfolio_eval = portfolio_eval,
      obj_value = objective_eval$value,
      obj_name = objective_eval$name,
      continuous_obj_value = continuous_objective_eval$value,
      pruned_obj_value = pruned_objective_eval$value,
      prune_obj_delta = prune_obj_delta,
      prune_obj_loss = prune_obj_loss,
      active_option_positions = active_option_positions,
      pruned_option_positions = pruned_option_positions,
      integer_obj_delta = integer_obj_delta,
      integer_obj_loss = integer_obj_loss,
      integer_candidate_count = integer_candidate_count,
      integer_feasible_count = integer_feasible_count,
      integer_used_initial_fallback = integer_used_initial_fallback,
      table_raw = table_result$raw,
      table_print = table_result$print,
      method_elapsed_seconds = method_elapsed_seconds
    )

    summary_row <- data.frame(
      constraint_mode = constraint_mode,
      objective = optimization_objective,
      utility_gamma = utility_gamma,
      obj_value = objective_eval$value,
      continuous_obj_value = continuous_objective_eval$value,
      pruned_obj_value = pruned_objective_eval$value,
      prune_obj_delta = prune_obj_delta,
      prune_obj_loss = prune_obj_loss,
      active_option_positions = active_option_positions,
      pruned_option_positions = pruned_option_positions,
      integer_obj_delta = integer_obj_delta,
      integer_obj_loss = integer_obj_loss,
      integer_candidate_count = integer_candidate_count,
      integer_feasible_count = integer_feasible_count,
      integer_used_initial_fallback = integer_used_initial_fallback,
      mean_variance_objective = portfolio_eval$mean_variance_objective,
      invested_weight = portfolio_eval$invested_weight,
      gross_position_weight = portfolio_eval$gross_position_weight,
      trade_cost = portfolio_eval$trade_cost,
      cash_after_trade = portfolio_eval$cash_after_trade,
      cash_fraction = portfolio_eval$cash_fraction,
      mean_wealth = portfolio_eval$mean_wealth,
      sd_wealth = portfolio_eval$sd_wealth,
      sharpe = portfolio_eval$sharpe,
      adj_sharpe = portfolio_eval$adj_sharpe,
      skew = portfolio_eval$skew,
      ex.kurt = portfolio_eval$ex.kurt,
      min_wealth = portfolio_eval$min_wealth,
      min_wealth_violation = portfolio_eval$min_wealth_violation,
      tail_slope = portfolio_eval$tail_slope,
      max_es_violation = portfolio_eval$max_es_violation,
      dependency_objective = portfolio_eval$dependency_objective,
      dependency_penalty = portfolio_eval$dependency_penalty,
      dependency_raw_min_eigen = dependency_diag$raw_min_eigen,
      dependency_raw_max_eigen = dependency_diag$raw_max_eigen,
      dependency_repaired_min_eigen = dependency_diag$repaired_min_eigen,
      dependency_repaired_max_eigen = dependency_diag$repaired_max_eigen,
      dependency_condition_number = dependency_diag$condition_number,
      dependency_psd_repaired = dependency_diag$repaired,
      portfolio_delta = portfolio_eval$portfolio_delta,
      portfolio_gamma = portfolio_eval$portfolio_gamma,
      portfolio_vega = portfolio_eval$portfolio_vega,
      delta_dollars = portfolio_eval$delta_dollars,
      gamma_dollars_1pct = portfolio_eval$gamma_dollars_1pct,
      vega_1pct = portfolio_eval$vega_1pct,
      delta_bound_violation = portfolio_eval$delta_bound_violation,
      current_mtm_wealth = portfolio_eval$current_mtm_wealth,
      mean_mtm_pnl = portfolio_eval$mean_mtm_pnl,
      sd_mtm_pnl = portfolio_eval$sd_mtm_pnl,
      max_mtm_var_violation = portfolio_eval$max_mtm_var_violation,
      max_mtm_es_violation = portfolio_eval$max_mtm_es_violation,
      hedged_mean_wealth = portfolio_eval$hedged_mean_wealth,
      hedged_sd_wealth = portfolio_eval$hedged_sd_wealth,
      hedged_sharpe = portfolio_eval$hedged_sharpe,
      hedged_skew = portfolio_eval$hedged_skew,
      hedged_ex.kurt = portfolio_eval$hedged_ex.kurt,
      hedged_mean_pnl = portfolio_eval$hedged_mean_pnl,
      hedged_sd_pnl = portfolio_eval$hedged_sd_pnl,
      hedged_avg_abs_final_delta = portfolio_eval$hedged_avg_abs_final_delta,
      optimizer = opt_diag$optimizer,
      optimizer_status = opt_diag$optimizer_status,
      optimizer_status_code = opt_diag$optimizer_status_code,
      optimizer_evals = opt_diag$optimizer_evals,
      optimizer_fallback_used = opt_diag$optimizer_fallback_used,
      method_elapsed_seconds = method_elapsed_seconds
    )
    for (i in seq_along(portfolio_eval$es_loss)) {
      summary_row[[portfolio_eval$es_names[i]]] <- portfolio_eval$es_loss[i]
    }
    for (i in seq_along(portfolio_eval$es_violation)) {
      summary_row[[portfolio_eval$es_violation_names[i]]] <- portfolio_eval$es_violation[i]
    }
    for (i in seq_along(portfolio_eval$mtm_var_values)) {
      summary_row[[portfolio_eval$mtm_var_names[i]]] <- portfolio_eval$mtm_var_values[i]
    }
    for (i in seq_along(portfolio_eval$mtm_var_violation)) {
      summary_row[[portfolio_eval$mtm_var_violation_names[i]]] <-
        portfolio_eval$mtm_var_violation[i]
    }
    for (i in seq_along(portfolio_eval$mtm_es_violation)) {
      summary_row[[portfolio_eval$mtm_es_violation_names[i]]] <-
        portfolio_eval$mtm_es_violation[i]
    }
    summary_rows <- rbind(summary_rows, summary_row)
    runs[[length(runs) + 1]] <- run
    if (print_runs) {
      print_option_optimization_run(run)
    }
  }
  }

  list(
    runs = runs,
    summary_rows = summary_rows,
    summary_print = format_option_optimization_summary(summary_rows)
  )
}

nonnegative_terminal_option_opt <- function(m,
                                            v,
                                            K,
                                            type,
                                            prices,
                                            expected_payoff,
                                            cov_payoff,
                                            dependency_matrix = NULL,
                                            budget,
                                            rf_growth,
                                            objective,
                                            risk_aversion,
                                            min_terminal_wealth = 0.0,
                                            initial_contracts = rep(0.0, length(K)),
                                            bid_prices = prices,
                                            ask_prices = prices,
                                            cvar_constraints = NULL,
                                            cvar_scenarios = NULL,
                                            payoff_scenarios = NULL,
                                            payoff_grid = NULL,
                                            lower_tail_slope = NULL,
                                            delta = NULL,
                                            portfolio_delta_bounds = c(NA_real_, NA_real_),
                                            max_abs_position_weight = 1.0,
                                            constrained_optimizer = "auto",
                                            max_starts = Inf,
                                            max_iter = 2000,
                                            extra_starts = list(),
                                            tol = 1e-10) {
  n <- length(K)
  if (is.null(payoff_grid)) {
    state_grid <- sort(unique(c(0.0, K)))
    payoff_grid <- option_payoff_matrix(state_grid, K, type)
  }
  tail_slope <- option_tail_slope(K, type)
  cvar_payoff_scenarios <- NULL
  if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) {
    cvar_payoff_scenarios <- option_payoff_matrix(cvar_scenarios, K, type)
  }
  cvar_tail_probs <- if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) cvar_constraints$tail_prob else numeric(0)
  cvar_max_loss <- if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) cvar_constraints$max_loss else numeric(0)
  cvar_payoff_scenarios_cpp <- if (is.null(cvar_payoff_scenarios)) {
    matrix(0.0, 0, n)
  } else {
    cvar_payoff_scenarios
  }
  delta_constraints <- portfolio_delta_constraint_rows(delta, portfolio_delta_bounds)

  ui <- rbind(
    payoff_grid - rf_growth * matrix(prices, nrow(payoff_grid), n, byrow = TRUE),
    tail_slope,
    if (!is.null(lower_tail_slope)) -lower_tail_slope else NULL,
    delta_constraints$ui,
    -prices,
    prices,
    -diag(prices),
    diag(prices)
  )
  ci <- c(
    rep(min_terminal_wealth - budget * rf_growth, nrow(payoff_grid)),
    0.0,
    if (!is.null(lower_tail_slope)) 0.0 else NULL,
    delta_constraints$ci,
    -budget,
    -budget,
    rep(-max_abs_position_weight * budget, n),
    rep(-max_abs_position_weight * budget, n)
  )

  objective_value <- function(contracts) {
    if (isTRUE(option_kernels_loaded) &&
        objective %in% c("mean_wealth", "mean_variance", "sharpe")) {
      return(nonnegative_terminal_objective_cpp(
        contracts = contracts,
        expected_payoff = expected_payoff,
        cov_payoff = cov_payoff,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        objective = objective,
        risk_aversion = risk_aversion,
        cvar_tail_probs = cvar_tail_probs,
        cvar_max_loss = cvar_max_loss,
        tol = tol
      ))
    }

    trade_cost <- rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
    excess <- sum(expected_payoff * contracts) - rf_growth * trade_cost
    variance <- as.numeric(t(contracts) %*% cov_payoff %*% contracts)
    sd <- sqrt(max(variance, 0.0))
    cvar_penalty <- 0.0
    if (!is.null(cvar_payoff_scenarios)) {
      cash <- budget - trade_cost
      terminal_wealth <- cash * rf_growth + as.numeric(cvar_payoff_scenarios %*% contracts)
      cvar_penalty <- cvar_constraint_penalty(terminal_wealth, budget, cvar_constraints)
    }

    if (objective == "mean_wealth") {
      return(excess - cvar_penalty)
    }

    if (objective == "mean_variance") {
      return(excess / budget - risk_aversion * variance / budget^2 - cvar_penalty)
    }

    if (objective == "dependency_penalty") {
      if (is.null(dependency_matrix)) {
        stop("dependency_penalty objective requires dependency_matrix")
      }
      position_weights <- position_exec_prices(contracts, bid_prices, ask_prices) *
        contracts / budget
      dependency_penalty <- as.numeric(t(position_weights) %*% dependency_matrix %*%
        position_weights)
      return(excess / budget - risk_aversion * dependency_penalty - cvar_penalty)
    }

    if (!is.finite(sd) || sd <= tol) {
      return(-Inf)
    }

    sharpe <- excess / sd

    if (objective == "sharpe") {
      return(sharpe - cvar_penalty)
    }

    if (objective == "adjusted_sharpe" || objective == "sharpe_adj") {
      stats <- if (is.null(payoff_scenarios)) {
        option_portfolio_stats(m, v, K, type, contracts)
      } else {
        option_portfolio_stats_from_scenarios(payoff_scenarios, contracts)
      }
      return(as.numeric(modified_sharpe_ratio(sharpe, stats["skew"], stats["ex.kurt"])) - cvar_penalty)
    }

    stop("unknown objective")
  }

  neg_objective <- function(contracts) {
    val <- objective_value(contracts)
    if (!is.finite(val)) {
      return(1e100)
    }
    -val
  }
  neg_objective_grad <- function(contracts) {
    grad <- nonnegative_terminal_gradient_cpp(
      contracts = contracts,
      expected_payoff = expected_payoff,
      cov_payoff = cov_payoff,
      cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
      initial_contracts = initial_contracts,
      bid_prices = bid_prices,
      ask_prices = ask_prices,
      budget = budget,
      rf_growth = rf_growth,
      objective = objective,
      risk_aversion = risk_aversion,
      cvar_tail_probs = cvar_tail_probs,
      cvar_max_loss = cvar_max_loss,
      tol = tol
    )
    if (any(!is.finite(grad))) {
      return(rep(0.0, length(contracts)))
    }
    -grad
  }
  constraint_values <- function(contracts) {
    if (isTRUE(option_kernels_loaded)) {
      return(option_constraint_values_cpp(
        contracts = contracts,
        ui = ui,
        ci = ci,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        cvar_tail_probs = cvar_tail_probs,
        cvar_max_loss = cvar_max_loss
      ))
    }

    linear_constraints <- as.numeric(ci - ui %*% contracts)
    if (is.null(cvar_payoff_scenarios)) {
      return(linear_constraints)
    }

    trade_cost <- rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
    cash <- budget - trade_cost
    terminal_wealth <- cash * rf_growth + as.numeric(cvar_payoff_scenarios %*% contracts)
    c(linear_constraints, cvar_constraint_excess(terminal_wealth, budget, cvar_constraints))
  }
  constraint_jacobian <- function(contracts) {
    if (isTRUE(option_kernels_loaded)) {
      return(option_constraint_jacobian_cpp(
        contracts = contracts,
        ui = ui,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        cvar_tail_probs = cvar_tail_probs
      ))
    }

    rbind(
      -ui,
      cvar_constraint_jacobian(
        contracts = contracts,
        budget = budget,
        rf_growth = rf_growth,
        cvar_constraints = cvar_constraints,
        cvar_payoff_scenarios = cvar_payoff_scenarios,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices
      )
    )
  }
  feasible_solution <- function(contracts, tol = 1e-6) {
    all(constraint_values(contracts) <= tol)
  }

  starts <- list(initial_contracts)
  first_call <- which(type == "call")[1]
  if (!is.na(first_call)) {
    x0 <- numeric(n)
    x0[first_call] <- 0.01 * budget / prices[first_call]
    starts[[length(starts) + 1]] <- x0
  }

  for (i in seq_len(n)) {
    x0 <- numeric(n)
    x0[i] <- 0.01 * budget / prices[i]
    starts <- add_contract_start(starts, x0, ui, ci)
  }
  for (x0 in extra_starts) {
    starts <- add_contract_start(starts, x0, ui, ci)
  }
  if (is.finite(max_starts) && length(starts) > max_starts) {
    starts <- starts[seq_len(max(1, max_starts))]
  }

  best_contracts <- starts[[1]]
  best_objective <- objective_value(best_contracts)
  best_optimizer <- "initial"
  best_status <- "initial"
  best_status_code <- NA_integer_
  best_evals <- 0L
  fallback_used <- FALSE

  use_nloptr <- constrained_optimizer %in% c("auto", "nloptr", "hybrid") &&
    isTRUE(option_kernels_loaded) &&
    objective %in% c("mean_wealth", "mean_variance", "sharpe") &&
    requireNamespace("nloptr", quietly = TRUE)
  use_nelder_mead <- constrained_optimizer %in% c("nelder_mead", "hybrid") ||
    !use_nloptr

  for (x0 in starts) {
    if (use_nloptr) {
      fit <- try(
        nloptr::nloptr(
          x0 = x0,
          eval_f = neg_objective,
          eval_grad_f = neg_objective_grad,
          eval_g_ineq = constraint_values,
          eval_jac_g_ineq = constraint_jacobian,
          opts = list(
            algorithm = "NLOPT_LD_SLSQP",
            maxeval = max_iter,
            xtol_rel = tol,
            print_level = 0
          )
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$solution)
        if (is.finite(val) && feasible_solution(fit$solution) && val > best_objective + 1e-12) {
          best_contracts <- fit$solution
          best_objective <- val
          best_optimizer <- "nloptr_slsqp"
          best_status <- decode_nloptr_status(fit$status)
          best_status_code <- as.integer(fit$status)
          best_evals <- as.integer(fit$iterations)
        }
      }
    }

    if (use_nelder_mead) {
      fit <- try(
        constrOptim(
          theta = x0,
          f = neg_objective,
          grad = NULL,
          ui = ui,
          ci = ci,
          method = "Nelder-Mead",
          control = list(maxit = max_iter)
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$par)
        if (is.finite(val) && val > best_objective + 1e-12) {
          best_contracts <- fit$par
          best_objective <- val
          best_optimizer <- "constrOptim_nelder_mead"
          best_status <- "ok"
          best_status_code <- NA_integer_
          best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
          fallback_used <- use_nloptr
        }
      }
    }
  }

  if (use_nloptr && !use_nelder_mead && best_optimizer == "initial") {
    fallback_used <- TRUE
    for (x0 in starts) {
      fit <- try(
        constrOptim(
          theta = x0,
          f = neg_objective,
          grad = NULL,
          ui = ui,
          ci = ci,
          method = "Nelder-Mead",
          control = list(maxit = max_iter)
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$par)
        if (is.finite(val) && val > best_objective + 1e-12) {
          best_contracts <- fit$par
          best_objective <- val
          best_optimizer <- "constrOptim_nelder_mead"
          best_status <- "fallback_ok"
          best_status_code <- NA_integer_
          best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
        }
      }
    }
  }

  list(
    contracts = best_contracts,
    weights = position_exec_prices(best_contracts, bid_prices, ask_prices) * best_contracts / budget,
    trade_contracts = rebalance_trade_contracts(best_contracts, initial_contracts),
    objective = best_objective,
    optimizer_diagnostics = list(
      optimizer = best_optimizer,
      optimizer_status = best_status,
      optimizer_status_code = best_status_code,
      optimizer_evals = best_evals,
      optimizer_fallback_used = fallback_used
    ),
    min_terminal_value = min(budget * rf_growth + as.numeric((payoff_grid - rf_growth *
      matrix(prices, nrow(payoff_grid), n, byrow = TRUE)) %*% best_contracts)),
    terminal_tail_slope = sum(tail_slope * best_contracts)
  )
}

expected_utility_value <- function(terminal_wealth, gamma) {
  if (any(terminal_wealth <= 0.0) || any(!is.finite(terminal_wealth))) {
    return(-Inf)
  }

  if (abs(gamma - 1.0) < 1e-10) {
    return(mean(log(terminal_wealth)))
  }

  mean(terminal_wealth^(1.0 - gamma) / (1.0 - gamma))
}

expected_utility_option_opt <- function(m,
                                        v,
                                        K,
                                        type,
                                        prices,
                                        budget,
                                        rf_growth,
                                        gamma,
                                        min_terminal_wealth,
                                        initial_contracts = rep(0.0, length(K)),
                                        bid_prices = prices,
                                        ask_prices = prices,
                                        cvar_constraints = NULL,
                                        cvar_scenarios = NULL,
                                        payoff_scenarios = NULL,
                                        payoff_grid = NULL,
                                        lower_tail_slope = NULL,
                                        delta = NULL,
                                        portfolio_delta_bounds = c(NA_real_, NA_real_),
                                        max_abs_position_weight = 1.0,
                                        max_starts = Inf,
                                        max_iter = 2000,
                                        constrained_optimizer = "auto",
                                        extra_starts = list(),
                                        tol = 1e-10,
                                        n_grid = 401) {
  n <- length(K)
  if (is.null(payoff_scenarios)) {
    z <- qnorm(seq(0.001, 0.999, length.out = n_grid))
    S <- exp(m + sqrt(v) * z)
    payoff_scenarios <- option_payoff_matrix(S, K, type)
  }
  cvar_payoff_scenarios <- NULL
  if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) {
    cvar_payoff_scenarios <- option_payoff_matrix(cvar_scenarios, K, type)
  }
  cvar_tail_probs <- if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) cvar_constraints$tail_prob else numeric(0)
  cvar_max_loss <- if (!is.null(cvar_constraints) && nrow(cvar_constraints) > 0) cvar_constraints$max_loss else numeric(0)
  cvar_payoff_scenarios_cpp <- if (is.null(cvar_payoff_scenarios)) {
    matrix(0.0, 0, n)
  } else {
    cvar_payoff_scenarios
  }
  if (is.null(payoff_grid)) {
    state_grid <- sort(unique(c(0.0, K)))
    payoff_grid <- option_payoff_matrix(state_grid, K, type)
  }
  tail_slope <- option_tail_slope(K, type)
  delta_constraints <- portfolio_delta_constraint_rows(delta, portfolio_delta_bounds)

  ui <- rbind(
    payoff_grid - rf_growth * matrix(prices, nrow(payoff_grid), n, byrow = TRUE),
    tail_slope,
    if (!is.null(lower_tail_slope)) -lower_tail_slope else NULL,
    delta_constraints$ui,
    -prices,
    prices,
    -diag(prices),
    diag(prices)
  )
  ci <- c(
    rep(min_terminal_wealth - budget * rf_growth, nrow(payoff_grid)),
    0.0,
    if (!is.null(lower_tail_slope)) 0.0 else NULL,
    delta_constraints$ci,
    -budget,
    -budget,
    rep(-max_abs_position_weight * budget, n),
    rep(-max_abs_position_weight * budget, n)
  )

  objective_value <- function(contracts) {
    if (isTRUE(option_kernels_loaded)) {
      return(expected_utility_objective_cpp(
        contracts = contracts,
        payoff_scenarios = payoff_scenarios,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        gamma = gamma,
        cvar_tail_probs = cvar_tail_probs,
        cvar_max_loss = cvar_max_loss
      ))
    }

    cash <- budget - rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
    terminal_wealth <- cash * rf_growth + as.numeric(payoff_scenarios %*% contracts)
    cvar_penalty <- 0.0
    if (!is.null(cvar_payoff_scenarios)) {
      cvar_wealth <- cash * rf_growth + as.numeric(cvar_payoff_scenarios %*% contracts)
      cvar_penalty <- cvar_constraint_penalty(cvar_wealth, budget, cvar_constraints)
    }
    expected_utility_value(terminal_wealth, gamma) - cvar_penalty
  }

  neg_objective <- function(contracts) {
    val <- objective_value(contracts)
    if (!is.finite(val)) {
      return(1e100)
    }
    -val
  }
  neg_objective_grad <- function(contracts) {
    grad <- expected_utility_gradient_cpp(
      contracts = contracts,
      payoff_scenarios = payoff_scenarios,
      cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
      initial_contracts = initial_contracts,
      bid_prices = bid_prices,
      ask_prices = ask_prices,
      budget = budget,
      rf_growth = rf_growth,
      gamma = gamma,
      cvar_tail_probs = cvar_tail_probs,
      cvar_max_loss = cvar_max_loss
    )
    if (any(!is.finite(grad))) {
      return(rep(0.0, length(contracts)))
    }
    -grad
  }
  constraint_values <- function(contracts) {
    if (isTRUE(option_kernels_loaded)) {
      return(option_constraint_values_cpp(
        contracts = contracts,
        ui = ui,
        ci = ci,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        cvar_tail_probs = cvar_tail_probs,
        cvar_max_loss = cvar_max_loss
      ))
    }

    linear_constraints <- as.numeric(ci - ui %*% contracts)
    if (is.null(cvar_payoff_scenarios)) {
      return(linear_constraints)
    }

    trade_cost <- rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
    cash <- budget - trade_cost
    terminal_wealth <- cash * rf_growth + as.numeric(cvar_payoff_scenarios %*% contracts)
    c(linear_constraints, cvar_constraint_excess(terminal_wealth, budget, cvar_constraints))
  }
  constraint_jacobian <- function(contracts) {
    if (isTRUE(option_kernels_loaded)) {
      return(option_constraint_jacobian_cpp(
        contracts = contracts,
        ui = ui,
        cvar_payoff_scenarios = cvar_payoff_scenarios_cpp,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        cvar_tail_probs = cvar_tail_probs
      ))
    }

    rbind(
      -ui,
      cvar_constraint_jacobian(
        contracts = contracts,
        budget = budget,
        rf_growth = rf_growth,
        cvar_constraints = cvar_constraints,
        cvar_payoff_scenarios = cvar_payoff_scenarios,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices
      )
    )
  }
  feasible_solution <- function(contracts, tol = 1e-6) {
    all(constraint_values(contracts) <= tol)
  }

  starts <- list(initial_contracts, numeric(n))
  first_call <- which(type == "call")[1]
  if (!is.na(first_call)) {
    x0 <- numeric(n)
    x0[first_call] <- min(0.01 * budget / prices[first_call], max_abs_position_weight * budget / prices[first_call])
    starts[[length(starts) + 1]] <- x0
  }

  for (i in seq_len(n)) {
    x0 <- numeric(n)
    x0[i] <- 0.01 * budget / prices[i]
    starts <- add_contract_start(starts, x0, ui, ci)
  }
  if (is.finite(max_starts) && length(starts) > max_starts) {
    starts <- starts[seq_len(max(1, max_starts))]
  }
  for (x0 in extra_starts) {
    starts <- add_contract_start(starts, x0, ui, ci)
  }

  best_contracts <- starts[[1]]
  best_objective <- objective_value(best_contracts)
  best_optimizer <- "initial"
  best_status <- "initial"
  best_status_code <- NA_integer_
  best_evals <- 0L
  fallback_used <- FALSE

  use_nloptr <- constrained_optimizer %in% c("auto", "nloptr", "hybrid") &&
    isTRUE(option_kernels_loaded) &&
    requireNamespace("nloptr", quietly = TRUE)
  use_nelder_mead <- constrained_optimizer %in% c("auto", "nelder_mead", "hybrid") &&
    !(use_nloptr && constrained_optimizer == "auto")

  for (x0 in starts) {
    if (use_nloptr) {
      fit <- try(
        nloptr::nloptr(
          x0 = x0,
          eval_f = neg_objective,
          eval_grad_f = neg_objective_grad,
          eval_g_ineq = constraint_values,
          eval_jac_g_ineq = constraint_jacobian,
          opts = list(
            algorithm = "NLOPT_LD_SLSQP",
            maxeval = max_iter,
            xtol_rel = tol,
            print_level = 0
          )
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$solution)
        if (is.finite(val) && feasible_solution(fit$solution) && val > best_objective + 1e-12) {
          best_contracts <- fit$solution
          best_objective <- val
          best_optimizer <- "nloptr_slsqp"
          best_status <- decode_nloptr_status(fit$status)
          best_status_code <- as.integer(fit$status)
          best_evals <- as.integer(fit$iterations)
        }
      }
    }

    if (use_nelder_mead) {
      fit <- try(
        constrOptim(
          theta = x0,
          f = neg_objective,
          grad = NULL,
          ui = ui,
          ci = ci,
          method = "Nelder-Mead",
          control = list(maxit = max_iter)
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$par)
        if (is.finite(val) && val > best_objective + 1e-12) {
          best_contracts <- fit$par
          best_objective <- val
          best_optimizer <- "constrOptim_nelder_mead"
          best_status <- "ok"
          best_status_code <- NA_integer_
          best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
          fallback_used <- use_nloptr
        }
      }
    }
  }

  if (use_nloptr && !use_nelder_mead && best_objective <= objective_value(starts[[1]]) + 1e-12) {
    fallback_used <- TRUE
    for (x0 in starts) {
      fit <- try(
        constrOptim(
          theta = x0,
          f = neg_objective,
          grad = NULL,
          ui = ui,
          ci = ci,
          method = "Nelder-Mead",
          control = list(maxit = max_iter)
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        val <- objective_value(fit$par)
        if (is.finite(val) && val > best_objective + 1e-12) {
          best_contracts <- fit$par
          best_objective <- val
          best_optimizer <- "constrOptim_nelder_mead"
          best_status <- "fallback_ok"
          best_status_code <- NA_integer_
          best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
        }
      }
    }
  }

  list(
    contracts = best_contracts,
    weights = position_exec_prices(best_contracts, bid_prices, ask_prices) * best_contracts / budget,
    trade_contracts = rebalance_trade_contracts(best_contracts, initial_contracts),
    objective = best_objective,
    optimizer_diagnostics = list(
      optimizer = best_optimizer,
      optimizer_status = best_status,
      optimizer_status_code = best_status_code,
      optimizer_evals = best_evals,
      optimizer_fallback_used = fallback_used
    )
  )
}

long_only_expected_utility_opt <- function(m,
                                           v,
                                           K,
                                           type,
                                           prices,
                                           budget,
                                           rf_growth,
                                           gamma,
                                           initial_contracts = rep(0.0, length(K)),
                                           bid_prices = prices,
                                           ask_prices = prices,
                                           payoff_scenarios = NULL,
                                           max_weight = 1.0,
                                           force_full_investment = FALSE,
                                           n_grid = 401) {
  n <- length(K)
  if (is.null(payoff_scenarios)) {
    z <- qnorm(seq(0.001, 0.999, length.out = n_grid))
    S <- exp(m + sqrt(v) * z)
    payoff_scenarios <- option_payoff_matrix(S, K, type)
  }

  softmax <- function(theta) {
    z <- theta - max(theta)
    ez <- exp(z)
    ez / sum(ez)
  }

  weights_from_theta <- function(theta) {
    if (force_full_investment) {
      return(max_weight * softmax(theta))
    }

    w_cash_and_options <- softmax(theta)
    max_weight * w_cash_and_options[-length(w_cash_and_options)]
  }

  objective_value <- function(weights) {
    contracts <- budget * weights / prices
    cash <- budget - rebalance_trade_cost(contracts, initial_contracts, bid_prices, ask_prices)
    terminal_wealth <- cash * rf_growth + as.numeric(payoff_scenarios %*% contracts)
    expected_utility_value(terminal_wealth, gamma)
  }

  objective <- function(theta) {
    val <- if (isTRUE(option_kernels_loaded)) {
      long_only_expected_utility_objective_cpp(
        theta = theta,
        payoff_scenarios = payoff_scenarios,
        prices = prices,
        initial_contracts = initial_contracts,
        bid_prices = bid_prices,
        ask_prices = ask_prices,
        budget = budget,
        rf_growth = rf_growth,
        gamma = gamma,
        max_weight = max_weight,
        force_full_investment = force_full_investment
      )
    } else {
      objective_value(weights_from_theta(theta))
    }
    if (!is.finite(val)) {
      return(1e100)
    }
    -val
  }
  objective_grad <- function(theta) {
    grad <- long_only_expected_utility_gradient_cpp(
      theta = theta,
      payoff_scenarios = payoff_scenarios,
      prices = prices,
      initial_contracts = initial_contracts,
      bid_prices = bid_prices,
      ask_prices = ask_prices,
      budget = budget,
      rf_growth = rf_growth,
      gamma = gamma,
      max_weight = max_weight,
      force_full_investment = force_full_investment
    )
    if (any(!is.finite(grad))) {
      return(rep(0.0, length(theta)))
    }
    -grad
  }

  theta_len <- if (force_full_investment) n else n + 1
  starts <- list(rep(0.0, theta_len))
  for (i in seq_len(theta_len)) {
    z0 <- rep(-8.0, theta_len)
    z0[i] <- 8.0
    starts[[length(starts) + 1]] <- z0
  }

  best_weights <- weights_from_theta(starts[[1]])
  best_objective <- objective_value(best_weights)
  best_optimizer <- "initial"
  best_status <- "initial"
  best_status_code <- NA_integer_
  best_evals <- 0L

  for (start in starts) {
    fit <- try(optim(
      start,
      objective,
      gr = if (isTRUE(option_kernels_loaded)) objective_grad else NULL,
      method = "BFGS",
      control = list(maxit = 1000)
    ), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      weights <- weights_from_theta(fit$par)
      val <- if (isTRUE(option_kernels_loaded)) {
        long_only_expected_utility_objective_cpp(
          theta = fit$par,
          payoff_scenarios = payoff_scenarios,
          prices = prices,
          initial_contracts = initial_contracts,
          bid_prices = bid_prices,
          ask_prices = ask_prices,
          budget = budget,
          rf_growth = rf_growth,
          gamma = gamma,
          max_weight = max_weight,
          force_full_investment = force_full_investment
        )
      } else {
        objective_value(weights)
      }
      if (is.finite(val) && val > best_objective + 1e-12) {
        best_weights <- weights
        best_objective <- val
        best_optimizer <- "optim_bfgs"
        best_status <- decode_optim_status(fit$convergence)
        best_status_code <- as.integer(fit$convergence)
        best_evals <- if (!is.null(fit$counts)) as.integer(sum(fit$counts)) else NA_integer_
      }
    }
  }

  list(
    weights = best_weights,
    contracts = budget * best_weights / prices,
    trade_contracts = rebalance_trade_contracts(
      budget * best_weights / prices,
      initial_contracts
    ),
    objective = best_objective,
    optimizer_diagnostics = list(
      optimizer = best_optimizer,
      optimizer_status = best_status,
      optimizer_status_code = best_status_code,
      optimizer_evals = best_evals,
      optimizer_fallback_used = FALSE
    )
  )
}

call_option_risk_table <- function(S0,
                                   m,
                                   v,
                                   K,
                                   implied_vol,
                                   r,
                                   q,
                                   T,
                                   weights) {
  stopifnot(length(K) == length(weights))
  if (length(implied_vol) == 1) {
    implied_vol <- rep(implied_vol, length(K))
  }

  stopifnot(length(K) == length(implied_vol))

  stats <- call_payoff_stats_table(m, v, K)
  prices <- bs_call_price_vec(S0, K, r, q, implied_vol, T)
  payoff_per_dollar <- stats$mean / prices
  sd_per_dollar <- stats$sd / prices
  edge <- payoff_per_dollar - exp(r * T)
  sharpe <- edge / sd_per_dollar
  adj_sharpe <- modified_sharpe_ratio(sharpe, stats$skew, stats$ex.kurt)

  data.frame(
    strike = K,
    implied_vol = implied_vol,
    price = as.numeric(prices),
    payoff_per_dollar = payoff_per_dollar,
    edge = edge,
    sharpe = sharpe,
    adj_sharpe = adj_sharpe,
    weight = weights,
    contracts = weights / prices
  )
}

call_portfolio_risk_stats <- function(m, v, K, contracts, initial_cost, r, T) {
  stats <- call_portfolio_stats(m, v, K, contracts)
  mean_payoff <- unname(stats["mean"])
  sd_payoff <- unname(stats["sd"])
  skew <- unname(stats["skew"])
  ex.kurt <- unname(stats["ex.kurt"])
  edge <- mean_payoff - initial_cost * exp(r * T)
  sharpe <- edge / sd_payoff
  adj_sharpe <- modified_sharpe_ratio(sharpe, skew, ex.kurt)

  c(
    initial_cost = initial_cost,
    payoff_per_dollar = mean_payoff / initial_cost,
    edge = edge / initial_cost,
    mean = mean_payoff,
    sd = sd_payoff,
    sharpe = sharpe,
    adj_sharpe = adj_sharpe,
    skew = skew,
    ex.kurt = ex.kurt
  )
}

# Raw moment of a weighted portfolio of call payoffs:
# sum_i weights[i] * (S_T - K[i])^+
call_portfolio_raw_moment <- function(m, v, K, weights, n) {
  stopifnot(length(K) == length(weights))

  stops <- sort(unique(K))
  total <- 0.0

  for (j in seq_along(stops)) {
    lower <- stops[j]
    if (j < length(stops)) {
      upper <- stops[j + 1]
    } else {
      upper <- Inf
    }

    active <- K <= lower
    a <- sum(weights[active])
    b <- sum(weights[active] * K[active])

    interval_value <- 0.0
    for (ell in 0:n) {
      interval_value <- interval_value +
        choose(n, ell) *
        a^ell *
        (-b)^(n - ell) *
        ln_interval_moment(m, v, lower, upper, ell)
    }

    total <- total + interval_value
  }

  total
}

# Mean, standard deviation, skewness, and excess kurtosis of a call portfolio
call_portfolio_stats <- function(m, v, K, weights) {
  r1 <- call_portfolio_raw_moment(m, v, K, weights, 1)
  r2 <- call_portfolio_raw_moment(m, v, K, weights, 2)
  r3 <- call_portfolio_raw_moment(m, v, K, weights, 3)
  r4 <- call_portfolio_raw_moment(m, v, K, weights, 4)

  mu2 <- r2 - r1^2
  mu3 <- r3 - 3 * r1 * r2 + 2 * r1^3
  mu4 <- r4 - 4 * r1 * r3 + 6 * r1^2 * r2 - 3 * r1^4

  c(
    mean = r1,
    sd = sqrt(mu2),
    skew = mu3 / mu2^(3 / 2),
    ex.kurt = mu4 / mu2^2 - 3
  )
}

# Table of mean, standard deviation, skewness, and excess kurtosis
call_payoff_stats_table <- function(m, v, K) {
  means <- sapply(K, function(k) call_payoff_mean(m, v, k))
  vars <- sapply(K, function(k) call_payoff_variance(m, v, k))
  sds <- sqrt(vars)
  skews <- sapply(K, function(k) call_payoff_skew(m, v, k))
  exkurts <- sapply(K, function(k) call_payoff_excess_kurtosis(m, v, k))

  out <- data.frame(
    strike = K,
    mean = means,
    sd = sds,
    skew = skews,
    ex.kurt = exkurts
  )

  rownames(out) <- NULL
  out
}

# Print summary for option payoff moments
print_call_payoff_summary <- function(m,
                                      v,
                                      K,
                                      print_variances = FALSE,
                                      print_covariances = FALSE,
                                      print_correlations = TRUE) {
  stats <- call_payoff_stats_table(m, v, K)
  covmat <- call_payoff_cov_mat(m, v, K)
  cormat <- call_payoff_cor_mat(covmat)

  cat("\nPayoff statistics:\n")
  print(round(stats, 3), row.names = FALSE)

  if (print_correlations) {
    cat("\nCorrelation matrix:\n")
    print(round(cormat, 3))
  }

  if (print_variances) {
    cat("\nVariances:\n")
    print(round(diag(covmat), 3))
  }

  if (print_covariances) {
    cat("\nCovariance matrix:\n")
    print(round(covmat, 3))
  }
}
