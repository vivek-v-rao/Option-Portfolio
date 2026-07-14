# Compare Corr(|X|, |Y|) for correlated bivariate normal and Student-t returns.

abs_corr_normal <- function(rho) {
  A <- (2 / pi) * (sqrt(1 - rho^2) + rho * asin(rho))
  return((A - 2 / pi) / (1 - 2 / pi))
}

abs_corr_student_t <- function(rho, nu) {
  if (nu <= 2) stop("nu must be > 2")
  A <- (2 / pi) * (sqrt(1 - rho^2) + rho * asin(rho))
  c_nu <- nu / (nu - 2)
  m_nu <- sqrt(2 / pi) * sqrt(nu / 2) *
    gamma((nu - 1) / 2) / gamma(nu / 2)
  return((c_nu * A - m_nu^2) / (c_nu - m_nu^2))
}

rho_grid <- seq(0.0, 0.8, by = 0.2)
nu_values <- c(30, 10, 5, 3)

out <- data.frame(
  rho = rho_grid,
  normal = abs_corr_normal(rho_grid)
)

for (nu in nu_values) {
  out[[paste0("t_df_", nu)]] <- abs_corr_student_t(rho_grid, nu)
}

print(round(out, 4), row.names = FALSE)
