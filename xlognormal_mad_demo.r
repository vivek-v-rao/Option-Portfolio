source("option_stats.r")

options(width = 10000)

lognormal_summary <- function(mu = 0.0, sigma_values = c(0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00)) {
  do.call(rbind, lapply(sigma_values, function(sigma) {
    v <- sigma^2
    mean_x <- exp(mu + 0.5 * v)
    variance <- (exp(v) - 1.0) * exp(2.0 * mu + v)
    skew <- (exp(v) + 2.0) * sqrt(exp(v) - 1.0)
    ex_kurt <- exp(4.0 * v) + 2.0 * exp(3.0 * v) + 3.0 * exp(2.0 * v) - 6.0

    data.frame(
      mu = mu,
      sigma = sigma,
      mean = mean_x,
      MAD = lognormal_mean_abs_deviation(mu, v),
      sd = sqrt(variance),
      "sd/MAD" = sqrt(variance) / lognormal_mean_abs_deviation(mu, v),
      skew = skew,
      ex_kurt = ex_kurt,
      check.names = FALSE
    )
  }))
}

format_lognormal_summary <- function(x) {
  x$mu <- sprintf("%.3f", x$mu)
  x$sigma <- sprintf("%.3f", x$sigma)
  x$mean <- sprintf("%.6f", x$mean)
  x$MAD <- sprintf("%.6f", x$MAD)
  x$sd <- sprintf("%.6f", x$sd)
  x[["sd/MAD"]] <- sprintf("%.6f", x[["sd/MAD"]])
  x$skew <- sprintf("%.6f", x$skew)
  x$ex_kurt <- sprintf("%.6f", x$ex_kurt)
  x
}

cat("MAD(x) = E(|x-mean(x)|)\n")
cat("Parameterization: log(X) ~ N(mu, sigma^2)\n\n")
cat("Normal distribution sd/MAD:", sprintf("%.6f", sqrt(pi / 2.0)), "\n\n")

cat("Lognormal distribution\n")
summary_table <- lognormal_summary()
print(format_lognormal_summary(summary_table), row.names = FALSE)
