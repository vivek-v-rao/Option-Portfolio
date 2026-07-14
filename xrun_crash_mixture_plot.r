option_scenario <- "crash_mixture"

option_config_overrides <- list(
	write_portfolio_plot = TRUE,
	portfolio_plot_file = "crash_mixture_plot.png",
	plot_constraint_mode = "long_only",
    plot_objective = "mean_variance"
)

source("xoptimize_options.r")
  