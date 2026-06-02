# precompute.R
# Run this script ONCE before launching app.R.
# It builds the full parameter grid, diagnoses all designs, and writes
# results.rds into shiny/.
#
# Usage (from the project root):
#   Rscript shiny/precompute.R
#
# Or from RStudio: open with the project root as working directory, then source.
# here::here() anchors all paths to the project root automatically.

source(here::here("shiny", "design_builder.R"))

set.seed(20260413)

# ============================================================================
# PARAMETER GRIDS - edit here to expand coverage
# ============================================================================

# Core grids shared by WP2 and WP3
N_grid <- c(2000, 4000, 6000)
# ## ADD: e.g. 1000, 3000, 5000 to fill gaps or extend range

tau_grid_wp2 <- c(0.05, 0.10, 0.20, 0.30)
# ## ADD: e.g. 0.08, 0.35, 0.40

# Keep WP3 survey-effect grid leaner by default for runtime.
tau_grid_wp3 <- c(0.05, 0.10, 0.20)
# ## ADD/RAISE: e.g. include 0.25, 0.30 for wider WP3 effect-size coverage

rho_grid <- c(0.20, 0.60)
# ## ADD: e.g. 0.10, 0.30, 0.50, 0.70

# WP2 attrition grids
attrition_grid <- c(0.10, 0.30, 0.50)
differential_attrition_grid <- c(0.00, 0.10)

# WP3 natural experiment grids
flood_exposure_grid <- c(0.10, 0.20)
# Share of respondents exposed to flood event between waves.

treat_prob_grid <- c(0.50, 0.65)
# Uneven assignment can increase precision for some estimands.

tau_flood_grid <- c(-0.10, -0.075, -0.05)
# Main effect of flood exposure on latent outcome.

tau_interaction_grid <- c(-0.08, 0.08)
# Interaction of survey experiment assignment and flood exposure (both signs).

flood_response_boost_grid <- c(0.00)
# Flood-exposure effect on follow-up retention, in log-odds units.
# Positive values increase follow-up response among flood-exposed; negative
# values reduce it.

# Keep WP3 panel attrition assumptions fixed by default for runtime,
# but keep as vectors so this is easy to customize later.
wp3_attrition_grid <- c(0.30)
wp3_differential_attrition_grid <- c(0.00)

# Simulations per design cell.
sims_precompute <- 200
# ## RAISE: to 500 or 1000 for final run

# -- Parallel setup -----------------------------------------------------------
use_parallel <- requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)

if (use_parallel) {
  library(future)
  library(future.apply)
  plan(multisession)
  n_cores <- future::availableCores() - 1
  message(sprintf(
    "[parallel] Running on %d cores via {future}. Expected speedup: ~%dx.",
    n_cores,
    n_cores
  ))
} else {
  message(
    "[sequential] {future} / {future.apply} not found - running single-threaded."
  )
  message("  To enable parallel execution:")
  message("  install.packages(c('future', 'future.apply'))")
  message("  Then re-run this script.")
}

# -- Build design grids -------------------------------------------------------

# WP2 post-only: attrition and natural-experiment parameters are irrelevant.
grid_wp2_post <- expand.grid(
  N = N_grid,
  tau = tau_grid_wp2,
  rho_y = rho_grid,
  stringsAsFactors = FALSE
)

# WP2 panel: vary standard panel attrition assumptions.
grid_wp2_panel <- expand.grid(
  N = N_grid,
  tau = tau_grid_wp2,
  rho_y = rho_grid,
  attrition_rate = attrition_grid,
  differential_attrition = differential_attrition_grid,
  stringsAsFactors = FALSE
)

# WP3 post-only: vary flood prevalence, treatment split, flood main effect,
# and interaction magnitude/sign.
grid_wp3_post <- expand.grid(
  N = N_grid,
  tau = tau_grid_wp3,
  rho_y = rho_grid,
  flood_exposure_rate = flood_exposure_grid,
  treat_prob = treat_prob_grid,
  tau_flood = tau_flood_grid,
  tau_interaction = tau_interaction_grid,
  stringsAsFactors = FALSE
)

# WP3 panel: same as WP3 post-only, plus flood-linked response boost and
# baseline panel attrition assumptions.
grid_wp3_panel <- expand.grid(
  N = N_grid,
  tau = tau_grid_wp3,
  rho_y = rho_grid,
  flood_exposure_rate = flood_exposure_grid,
  treat_prob = treat_prob_grid,
  tau_flood = tau_flood_grid,
  tau_interaction = tau_interaction_grid,
  attrition_rate = wp3_attrition_grid,
  differential_attrition = wp3_differential_attrition_grid,
  flood_response_boost = flood_response_boost_grid,
  stringsAsFactors = FALSE
)

total_cells <- nrow(grid_wp2_post) +
  nrow(grid_wp2_panel) +
  nrow(grid_wp3_post) +
  nrow(grid_wp3_panel)

message(sprintf(
  paste0(
    "\nGrid summary",
    "\n  WP2 post-only : %d cells",
    "\n  WP2 panel     : %d cells",
    "\n  WP3 post-only : %d cells",
    "\n  WP3 panel     : %d cells",
    "\n  Total         : %d cells x %d sims = %s simulations\n"
  ),
  nrow(grid_wp2_post),
  nrow(grid_wp2_panel),
  nrow(grid_wp3_post),
  nrow(grid_wp3_panel),
  total_cells,
  sims_precompute,
  format(total_cells * sims_precompute, big.mark = ",")
))

# -- Build design objects -----------------------------------------------------
message("Building WP2 post-only designs...")
wp2_post_designs <- lapply(seq_len(nrow(grid_wp2_post)), function(i) {
  build_wp2_design(
    N = grid_wp2_post$N[i],
    tau = grid_wp2_post$tau[i],
    rho_y = grid_wp2_post$rho_y[i],
    attrition_rate = 0.30,
    differential_attrition = 0.00,
    flood_exposure_rate = 0.10,
    treat_prob = 0.50,
    tau_flood = -0.10,
    tau_interaction = 0.00,
    flood_response_boost = 0.00,
    type = "wp2_post_only"
  )
})
names(wp2_post_designs) <- paste0("wp2po_", seq_len(nrow(grid_wp2_post)))

message("Building WP2 panel designs...")
wp2_panel_designs <- lapply(seq_len(nrow(grid_wp2_panel)), function(i) {
  build_wp2_design(
    N = grid_wp2_panel$N[i],
    tau = grid_wp2_panel$tau[i],
    rho_y = grid_wp2_panel$rho_y[i],
    attrition_rate = grid_wp2_panel$attrition_rate[i],
    differential_attrition = grid_wp2_panel$differential_attrition[i],
    flood_exposure_rate = 0.10,
    treat_prob = 0.50,
    tau_flood = -0.10,
    tau_interaction = 0.00,
    flood_response_boost = 0.00,
    type = "wp2_panel"
  )
})
names(wp2_panel_designs) <- paste0("wp2panel_", seq_len(nrow(grid_wp2_panel)))

message("Building WP3 post-only designs...")
wp3_post_designs <- lapply(seq_len(nrow(grid_wp3_post)), function(i) {
  build_wp2_design(
    N = grid_wp3_post$N[i],
    tau = grid_wp3_post$tau[i],
    rho_y = grid_wp3_post$rho_y[i],
    flood_exposure_rate = grid_wp3_post$flood_exposure_rate[i],
    treat_prob = grid_wp3_post$treat_prob[i],
    tau_flood = grid_wp3_post$tau_flood[i],
    tau_interaction = grid_wp3_post$tau_interaction[i],
    attrition_rate = 0.30,
    differential_attrition = 0.00,
    flood_response_boost = 0.00,
    type = "wp3_post_only"
  )
})
names(wp3_post_designs) <- paste0("wp3po_", seq_len(nrow(grid_wp3_post)))

message("Building WP3 panel designs...")
wp3_panel_designs <- lapply(seq_len(nrow(grid_wp3_panel)), function(i) {
  build_wp2_design(
    N = grid_wp3_panel$N[i],
    tau = grid_wp3_panel$tau[i],
    rho_y = grid_wp3_panel$rho_y[i],
    flood_exposure_rate = grid_wp3_panel$flood_exposure_rate[i],
    treat_prob = grid_wp3_panel$treat_prob[i],
    tau_flood = grid_wp3_panel$tau_flood[i],
    tau_interaction = grid_wp3_panel$tau_interaction[i],
    attrition_rate = grid_wp3_panel$attrition_rate[i],
    differential_attrition = grid_wp3_panel$differential_attrition[i],
    flood_response_boost = grid_wp3_panel$flood_response_boost[i],
    type = "wp3_panel"
  )
})
names(wp3_panel_designs) <- paste0("wp3panel_", seq_len(nrow(grid_wp3_panel)))

# -- Diagnose -----------------------------------------------------------------
t0 <- proc.time()

message(sprintf(
  "\nDiagnosing WP2 post-only (%d cells)...",
  nrow(grid_wp2_post)
))
diag_wp2_post <- diagnose_designs(
  wp2_post_designs,
  diagnosands = diagnosands_power,
  sims = sims_precompute,
  bootstrap_sims = FALSE
)

message(sprintf("Diagnosing WP2 panel (%d cells)...", nrow(grid_wp2_panel)))
diag_wp2_panel <- diagnose_designs(
  wp2_panel_designs,
  diagnosands = diagnosands_power,
  sims = sims_precompute,
  bootstrap_sims = FALSE
)

message(sprintf("Diagnosing WP3 post-only (%d cells)...", nrow(grid_wp3_post)))
diag_wp3_post <- diagnose_designs(
  wp3_post_designs,
  diagnosands = diagnosands_power,
  sims = sims_precompute,
  bootstrap_sims = FALSE
)

message(sprintf("Diagnosing WP3 panel (%d cells)...", nrow(grid_wp3_panel)))
diag_wp3_panel <- diagnose_designs(
  wp3_panel_designs,
  diagnosands = diagnosands_power,
  sims = sims_precompute,
  bootstrap_sims = FALSE
)

elapsed_min <- (proc.time() - t0)["elapsed"] / 60
message(sprintf("\nDiagnosis complete in %.1f minutes.", elapsed_min))

# -- Join parameter values back to diagnosands --------------------------------
parse_idx <- function(design_col, prefix) {
  as.integer(sub(prefix, "", as.character(design_col), fixed = TRUE))
}

results_wp2_post <- diag_wp2_post$diagnosands_df |>
  mutate(
    design_idx = parse_idx(design, "wp2po_"),
    design_type = "wp2_post_only"
  ) |>
  left_join(
    grid_wp2_post |>
      mutate(
        design_idx = row_number(),
        attrition_rate = NA_real_,
        differential_attrition = NA_real_,
        flood_exposure_rate = NA_real_,
        treat_prob = NA_real_,
        tau_flood = NA_real_,
        tau_interaction = NA_real_,
        flood_response_boost = NA_real_
      ),
    by = "design_idx"
  )

results_wp2_panel <- diag_wp2_panel$diagnosands_df |>
  mutate(
    design_idx = parse_idx(design, "wp2panel_"),
    design_type = "wp2_panel"
  ) |>
  left_join(
    grid_wp2_panel |>
      mutate(
        design_idx = row_number(),
        flood_exposure_rate = NA_real_,
        treat_prob = NA_real_,
        tau_flood = NA_real_,
        tau_interaction = NA_real_,
        flood_response_boost = NA_real_
      ),
    by = "design_idx"
  )

results_wp3_post <- diag_wp3_post$diagnosands_df |>
  mutate(
    design_idx = parse_idx(design, "wp3po_"),
    design_type = "wp3_post_only"
  ) |>
  left_join(
    grid_wp3_post |>
      mutate(
        design_idx = row_number(),
        attrition_rate = NA_real_,
        differential_attrition = NA_real_,
        flood_response_boost = NA_real_
      ),
    by = "design_idx"
  )

results_wp3_panel <- diag_wp3_panel$diagnosands_df |>
  mutate(
    design_idx = parse_idx(design, "wp3panel_"),
    design_type = "wp3_panel"
  ) |>
  left_join(
    grid_wp3_panel |>
      mutate(design_idx = row_number()),
    by = "design_idx"
  )

label_lookup <- setNames(
  vapply(DESIGN_REGISTRY, function(x) x$label, character(1)),
  vapply(DESIGN_REGISTRY, function(x) x$type, character(1))
)

results <- bind_rows(
  results_wp2_post,
  results_wp2_panel,
  results_wp3_post,
  results_wp3_panel
) |>
  select(-design_idx) |>
  mutate(
    estimator_label = dplyr::recode(estimator, !!!ESTIMATOR_LABELS),
    outcome_scale = ifelse(
      grepl("1_5", estimator),
      "Likert 1-5",
      "Latent SD"
    ),
    design_label_display = dplyr::recode(design_type, !!!label_lookup)
  )

# -- Save ---------------------------------------------------------------------
# Store grid metadata as an attribute so the app knows which values are valid.
attr(results, "grid_meta") <- list(
  N_grid = N_grid,
  tau_grid_wp2 = tau_grid_wp2,
  tau_grid_wp3 = tau_grid_wp3,
  tau_grid = sort(unique(c(tau_grid_wp2, tau_grid_wp3))),
  rho_grid = rho_grid,
  attrition_grid = attrition_grid,
  differential_attrition_grid = differential_attrition_grid,
  flood_exposure_grid = flood_exposure_grid,
  treat_prob_grid = treat_prob_grid,
  tau_flood_grid = tau_flood_grid,
  tau_interaction_grid = tau_interaction_grid,
  flood_response_boost_grid = flood_response_boost_grid,
  wp3_attrition_grid = wp3_attrition_grid,
  wp3_differential_attrition_grid = wp3_differential_attrition_grid,
  sims = sims_precompute,
  timestamp = Sys.time()
)

out_path <- here::here("shiny", "results.rds")
saveRDS(results, out_path)
message(sprintf(
  "\nSaved: %s  (%d rows, %d columns)",
  out_path,
  nrow(results),
  ncol(results)
))
message("Grid metadata: attr(readRDS(out_path), 'grid_meta')")
message("\nNext step: shiny::runApp(here::here('shiny', 'app.R'))")
