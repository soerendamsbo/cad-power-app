# findings_grid.R
#
# Comprehensive parameter sweep for updating "Preliminary Findings".
# Covers WP3 post-only and WP3 panel designs with expanded grids designed
# to answer four key questions:
#   1. At what (N, flood_exposure_rate) combinations is flood-effect power
#      adequate?
#   2. Under what (rho_y, attrition) combinations does panel clearly beat
#      post-only?
#   3. Is interaction power ever adequate, and for what (tau_interaction, N)?
#   4. How sensitive is panel power to flood-linked retention bias
#      (flood_response_boost)?
#
# Usage (from project root):
#   Rscript shiny/findings_grid.R
#   Rscript shiny/findings_grid.R sims=200 out_dir=shiny/findings_grid_out
#
# Outputs (written to out_dir):
#   findings_grid.csv  — one row per (cell × estimator × inquiry), all
#                        parameters + diagnosands
#   findings_grid.rds  — list(results = <data.frame>, meta = <list>)

suppressPackageStartupMessages({
  library(dplyr)
})

# -- Argument parsing ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  key <- paste0(name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(key, "", hit[[1]], fixed = TRUE)
}

sims <- as.integer(get_arg("sims", "200"))
out_dir_arg <- get_arg("out_dir", "shiny/findings_grid_out")

if (is.na(sims) || sims < 20) {
  stop("`sims` must be an integer >= 20.")
}

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package `here` is required. Install with install.packages('here').")
}

source(here::here("shiny", "design_builder.R"))

out_dir <- if (grepl("^/", out_dir_arg)) {
  out_dir_arg
} else {
  here::here(out_dir_arg)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260413)

# -- Parallel setup ------------------------------------------------------------
use_parallel <- requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)

if (use_parallel) {
  library(future)
  library(future.apply)
  plan(multisession)
  n_cores <- future::availableCores() - 1
  message(sprintf("[parallel] %d cores via {future}.", n_cores))
} else {
  message(
    "[sequential] Install {future} + {future.apply} for parallel execution."
  )
}

# -- Parameter grids -----------------------------------------------------------

# WP3 post-only: rho_y is irrelevant (unit-variance outcome regardless of rho_y
# when not using baseline covariate). Fixed at 0.40 for model consistency.
# treat_prob fixed at 0.50; prior quick sweep showed 0.50 vs 0.65 had minimal
# impact on flood-effect power, and unequal splits are operationally constrained.
grid_post <- expand.grid(
  N = c(2000, 3000, 4000, 6000),
  tau = c(0.05, 0.10, 0.20),
  flood_exposure_rate = c(0.05, 0.10, 0.15, 0.20),
  tau_flood = c(-0.15, -0.10, -0.05),
  tau_interaction = c(-0.08, 0.00, 0.08),
  stringsAsFactors = FALSE
)

# WP3 panel: adds rho_y (drives baseline-adjustment gain), attrition_rate
# (erodes effective N), and flood_response_boost (composition sensitivity).
grid_panel <- expand.grid(
  N = c(2000, 4000, 6000),
  tau = c(0.05, 0.10, 0.20),
  rho_y = c(0.20, 0.40, 0.60),
  attrition_rate = c(0.10, 0.30, 0.50),
  flood_exposure_rate = c(0.05, 0.10, 0.20),
  tau_flood = c(-0.15, -0.10, -0.05),
  tau_interaction = c(-0.08, 0.00, 0.08),
  flood_response_boost = c(0.00, 0.05),
  stringsAsFactors = FALSE
)

total_cells <- nrow(grid_post) + nrow(grid_panel)
message(sprintf(
  paste0(
    "\nGrid summary",
    "\n  WP3 post-only : %d cells",
    "\n  WP3 panel     : %d cells",
    "\n  Total         : %d cells x %d sims = %s simulations"
  ),
  nrow(grid_post),
  nrow(grid_panel),
  total_cells,
  sims,
  format(total_cells * sims, big.mark = ",")
))

# -- Helpers: run one design cell and attach its parameter values --------------

run_post <- function(i) {
  g <- grid_post[i, ]
  d <- build_wp2_design(
    N = g$N,
    tau = g$tau,
    tau_flood = g$tau_flood,
    tau_interaction = g$tau_interaction,
    rho_y = 0.40,
    attrition_rate = 0.30,
    differential_attrition = 0.00,
    flood_exposure_rate = g$flood_exposure_rate,
    treat_prob = 0.50,
    flood_response_boost = 0.00,
    type = "wp3_post_only"
  )
  diagnose_design(
    d,
    diagnosands = diagnosands_power,
    sims = sims,
    bootstrap_sims = FALSE
  )$diagnosands_df |>
    mutate(
      design_type = "wp3_post_only",
      N = g$N,
      tau = g$tau,
      rho_y = NA_real_,
      attrition_rate = NA_real_,
      flood_exposure_rate = g$flood_exposure_rate,
      tau_flood = g$tau_flood,
      tau_interaction = g$tau_interaction,
      flood_response_boost = NA_real_,
      treat_prob = 0.50
    )
}

run_panel <- function(i) {
  g <- grid_panel[i, ]
  d <- build_wp2_design(
    N = g$N,
    tau = g$tau,
    tau_flood = g$tau_flood,
    tau_interaction = g$tau_interaction,
    rho_y = g$rho_y,
    attrition_rate = g$attrition_rate,
    differential_attrition = 0.00,
    flood_exposure_rate = g$flood_exposure_rate,
    treat_prob = 0.50,
    flood_response_boost = g$flood_response_boost,
    type = "wp3_panel"
  )
  diagnose_design(
    d,
    diagnosands = diagnosands_power,
    sims = sims,
    bootstrap_sims = FALSE
  )$diagnosands_df |>
    mutate(
      design_type = "wp3_panel",
      N = g$N,
      tau = g$tau,
      rho_y = g$rho_y,
      attrition_rate = g$attrition_rate,
      flood_exposure_rate = g$flood_exposure_rate,
      tau_flood = g$tau_flood,
      tau_interaction = g$tau_interaction,
      flood_response_boost = g$flood_response_boost,
      treat_prob = 0.50
    )
}

# -- Run diagnoses -------------------------------------------------------------
t0 <- proc.time()

message(sprintf("\nRunning WP3 post-only (%d cells)...", nrow(grid_post)))
if (use_parallel) {
  post_res <- bind_rows(
    future_lapply(seq_len(nrow(grid_post)), run_post, future.seed = TRUE)
  )
} else {
  post_res <- bind_rows(lapply(seq_len(nrow(grid_post)), run_post))
}

message(sprintf("Running WP3 panel (%d cells)...", nrow(grid_panel)))
if (use_parallel) {
  panel_res <- bind_rows(
    future_lapply(seq_len(nrow(grid_panel)), run_panel, future.seed = TRUE)
  )
} else {
  panel_res <- bind_rows(lapply(seq_len(nrow(grid_panel)), run_panel))
}

elapsed_min <- (proc.time() - t0)["elapsed"] / 60
message(sprintf("\nDiagnosis complete in %.1f minutes.", elapsed_min))

# -- Combine and annotate ------------------------------------------------------
results <- bind_rows(post_res, panel_res) |>
  mutate(
    estimator_label = dplyr::recode(estimator, !!!ESTIMATOR_LABELS),
    outcome_scale = ifelse(grepl("1_5", estimator), "Likert 1-5", "Latent SD"),
    expected_flooded_n = N * flood_exposure_rate
  )

# -- Save outputs --------------------------------------------------------------
run_meta <- list(
  timestamp = Sys.time(),
  sims = sims,
  seed = 20260413,
  total_cells = total_cells,
  total_sims = total_cells * sims,
  elapsed_minutes = elapsed_min,
  grids = list(
    post = grid_post,
    panel = grid_panel
  )
)

out <- list(results = results, meta = run_meta)
saveRDS(out, file.path(out_dir, "findings_grid.rds"))
write.csv(results, file.path(out_dir, "findings_grid.csv"), row.names = FALSE)

message(sprintf(
  "\nSaved to: %s\n  findings_grid.rds\n  findings_grid.csv  (%d rows, %d columns)",
  out_dir,
  nrow(results),
  ncol(results)
))

message("\nColumn reference:")
message("  design_type, N, tau, rho_y, attrition_rate, flood_exposure_rate,")
message("  tau_flood, tau_interaction, flood_response_boost, treat_prob,")
message("  estimator, inquiry, term, outcome_scale,")
message("  power, bias, rmse, coverage, mean_estimate, mean_estimand,")
message("  mean_se, type_s_error, n_sims, expected_flooded_n")

message("\nKey inquiry values for filtering:")
message("  flood main effect  : inquiry == 'flood_avg_latent'")
message("  flood interaction  : inquiry == 'flood_survey_contrast_latent'")
message("  survey main effect : inquiry == 'survey_avg_latent'")
message(
  "  subgroups          : inquiry %in% c('flood_when_survey_control_latent',"
)
message("                         'flood_when_survey_treated_latent')")
