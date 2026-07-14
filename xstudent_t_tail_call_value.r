# Study call prices when Student-t log returns are modified to have a finite
# right tail. Each modified distribution is recentered so E[S_T] equals the
# forward price S0 * exp((r - q) * T).

source("option_stats.r")
options(width = 10000)

S0 <- 100
K_over_S0 <- c(1.0, 1.2, 1.4)
K_values <- S0 * K_over_S0
r <- 0.03
q <- 0.00
T <- 1.0
realized_sigma <- 0.20
df <- 5.0
cutoffs <- seq(0.2, 1.2, by = 0.1)
exp_tail_lambdas <- c(1.25, 1.5, 2.0)
gaussian_tail_lambdas <- c(0.5, 1.0, 2.0)
csv_file <- "outputs/student_t_tail_call_prices.csv"

cat("Student-t right-tail treatment call study\n")
cat("S0:", S0, "\n")
cat("K_over_S0:", K_over_S0, "\n")
cat("K_values:", K_values, "\n")
cat("r:", r, "\n")
cat("q:", q, "\n")
cat("T:", T, "\n")
cat("realized_sigma:", realized_sigma, "\n")
cat("df:", df, "\n")
cat("cutoffs:", cutoffs, "\n")
cat("exp_tail_lambdas:", exp_tail_lambdas, "\n")
cat("gaussian_tail_lambdas:", gaussian_tail_lambdas, "\n")
black_scholes_prices <- data.frame(
  K_over_S0 = K_over_S0,
  K = K_values,
  black_scholes_call_price = bs_call_price_vec(S0, K_values, r, q, realized_sigma, T)
)
black_scholes_prices$black_scholes_call_price <- round(
  black_scholes_prices$black_scholes_call_price,
  6
)
cat("Black-Scholes call prices\n")
print(black_scholes_prices, row.names = FALSE)
cat("\n")

results <- do.call(rbind, lapply(seq_along(K_values), function(i) {
  one <- student_t_tail_call_price_table(
    S0 = S0, K = K_values[i], r = r, q = q, T = T,
    realized_sigma = realized_sigma, df = df, cutoffs = cutoffs,
    exp_tail_lambdas = exp_tail_lambdas,
    gaussian_tail_lambdas = gaussian_tail_lambdas
  )
  one$K_over_S0 <- K_over_S0[i]
  one$K <- K_values[i]
  one$black_scholes_call_price <- black_scholes_prices$black_scholes_call_price[i]
  one
}))

results$method <- results$treatment
results$method[results$treatment == "exp_tail"] <- paste0(
  "exp_tail_", results$exp_tail_lambda[results$treatment == "exp_tail"]
)
results$method[results$treatment == "gaussian_tail"] <- paste0(
  "gaussian_tail_", results$gaussian_tail_lambda[results$treatment == "gaussian_tail"]
)

price_table <- reshape(
  results[, c("K_over_S0", "K", "cutoff", "method", "call_price")],
  idvar = c("K_over_S0", "K", "cutoff"),
  timevar = "method",
  direction = "wide"
)
names(price_table) <- sub("^call_price\\.", "", names(price_table))
price_table <- price_table[order(price_table$K_over_S0, price_table$cutoff), ]

format_price_table <- price_table
numeric_cols <- names(format_price_table)[sapply(format_price_table, is.numeric)]
format_price_table[numeric_cols] <- lapply(format_price_table[numeric_cols], function(x) round(x, 4))

cat("Call price by strike, cutoff, and tail treatment\n")
print(format_price_table, row.names = FALSE)
cat("\n")

diagnostics <- results[
  results$cutoff %in% c(min(cutoffs), median(cutoffs), max(cutoffs)),
  c(
    "K_over_S0", "K", "cutoff", "method", "recenter_shift",
    "exp_moment_before_shift", "normalization", "atom_mass",
    "tail_prob_after_recenter", "forward_check", "call_price",
    "black_scholes_call_price"
  )
]
diagnostics[, c(
  "recenter_shift", "exp_moment_before_shift", "normalization", "atom_mass",
  "tail_prob_after_recenter", "forward_check", "call_price",
  "black_scholes_call_price"
)] <- lapply(
  diagnostics[, c(
    "recenter_shift", "exp_moment_before_shift", "normalization", "atom_mass",
    "tail_prob_after_recenter", "forward_check", "call_price",
    "black_scholes_call_price"
  )],
  function(x) round(x, 6)
)

cat("Selected diagnostics\n")
print(diagnostics, row.names = FALSE)
cat("\n")

dir <- dirname(csv_file)
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
write.csv(results, csv_file, row.names = FALSE)
cat("csv_file:", csv_file, "\n")
