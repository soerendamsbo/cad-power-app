# quick_findings_sweep.R
#
# Purpose:
#   Recompute the lightweight WP3 "preliminary findings" snapshot used for
#   first-impression planning in the app.
#
# Usage (from project root):
#   Rscript shiny/quick_findings_sweep.R
#   Rscript shiny/quick_findings_sweep.R sims=80
#   Rscript shiny/quick_findings_sweep.R sims=80 out_dir=shiny/quick_findings
#
# Notes:
# - This is intentionally a focused, low-cost sweep for directional guidance.
# - It does not replace full precompute runs.

suppressPackageStartupMessages({
  library(dplyr)
})

# -- argument parsing ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  key <- paste0(name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[1]], fixed = TRUE)
}

sims <- as.integer(get_arg("sims", "40"))
out_dir_arg <- get_arg("out_dir", "shiny/quick_findings")

if (is.na(sims) || sims < 20) {
  stop("`sims` must be an integer >= 20.")
}

# -- locate and source design builder -----------------------------------------
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package `here` is required. Install with install.packages('here').")
}

builder_path <- here::here("shiny", "design_builder.R")
if (!file.exists(builder_path)) {
  stop("Could not find design_builder.R at: ", builder_path)
}
source(builder_path)

out_dir <- if (grepl("^/", out_dir_arg)) out_dir_arg else here::here(out_dir_arg)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- helper -------------------------------------------------------------------
run_one <- function(
  type,
  N,
  rho_y,
  attrition_rate,
  treat_prob,
  flood_exposure_rate,
  tau_flood,
  tau_interaction,
  sims_local
) {
  d <- build_wp2_design(
    N = N,
    tau = 0.10,
    tau_flood = tau_flood,
    tau_interaction = tau_interaction,
    rho_y = rho_y,
    attrition_rate = attrition_rate,
    differential_attrition = 0.00,
    flood_exposure_rate = flood_exposure_rate,
    treat_prob = treat_prob,
    flood_response_boost = 0.05,
    type = type
  )

  diagnose_design(
    d,
    diagnosands = diagnosands_power,
    sims = sims_local,
    bootstrap_sims = FALSE
  )$diagnosands_df |>
    mutate(
      N = N,
      rho_y = rho_y,
      attrition_rate = attrition_rate,
      treat_prob = treat_prob,
      flood_exposure_rate = flood_exposure_rate,
      tau_flood = tau_flood,
      tau_interaction = tau_interaction,
      design_type = type
    )
}

# -- focused sweep grid --------------------------------------------------------
Ns <- c(2000, 4000)
rhos <- c(0.20, 0.60)
attrs <- c(0.10, 0.50)
treats <- c(0.50, 0.65)
floods <- c(0.10, 0.20)
tau_floods <- c(-0.10, -0.05)
tau_inters <- c(-0.04, 0.04)

post_grid <- expand.grid(
  N = Ns,
  rho_y = rhos,
  treat_prob = treats,
  flood_exposure_rate = floods,
  tau_flood = tau_floods,
  tau_interaction = tau_inters,
  stringsAsFactors = FALSE
)

panel_grid <- expand.grid(
  N = Ns,
  rho_y = rhos,
  attrition_rate = attrs,
  treat_prob = treats,
  flood_exposure_rate = floods,
  tau_flood = tau_floods,
  tau_interaction = tau_inters,
  stringsAsFactors = FALSE
)

message(sprintf("Running WP3 post-only cells: %d", nrow(post_grid)))
post_res <- bind_rows(lapply(seq_len(nrow(post_grid)), function(i) {
  g <- post_grid[i, ]
  run_one(
    type = "wp3_post_only",
    N = g$N,
    rho_y = g$rho_y,
    attrition_rate = 0.30,
    treat_prob = g$treat_prob,
    flood_exposure_rate = g$flood_exposure_rate,
    tau_flood = g$tau_flood,
    tau_interaction = g$tau_interaction,
    sims_local = sims
  )
}))

message(sprintf("Running WP3 panel cells: %d", nrow(panel_grid)))
panel_res <- bind_rows(lapply(seq_len(nrow(panel_grid)), function(i) {
  g <- panel_grid[i, ]
  run_one(
    type = "wp3_panel",
    N = g$N,
    rho_y = g$rho_y,
    attrition_rate = g$attrition_rate,
    treat_prob = g$treat_prob,
    flood_exposure_rate = g$flood_exposure_rate,
    tau_flood = g$tau_flood,
    tau_interaction = g$tau_interaction,
    sims_local = sims
  )
}))

all_res <- bind_rows(post_res, panel_res)

# -- summaries used by Preliminary findings -----------------------------------
flood_power <- all_res |>
  filter(estimator %in% c("wp3_post_flood_latent", "wp3_panel_flood_latent")) |>
  mutate(panel = grepl("panel", estimator))

inter_power <- all_res |>
  filter(estimator %in% c("wp3_post_interaction_latent", "wp3_panel_interaction_latent")) |>
  mutate(panel = grepl("panel", estimator))

post_cmp <- flood_power |>
  filter(!panel) |>
  transmute(
    N,
    rho_y,
    treat_prob,
    flood_exposure_rate,
    tau_flood,
    tau_interaction,
    post_power = power
  )

panel_cmp <- flood_power |>
  filter(panel) |>
  transmute(
    N,
    rho_y,
    attrition_rate,
    treat_prob,
    flood_exposure_rate,
    tau_flood,
    tau_interaction,
    panel_power = power
  )

joined <- inner_join(
  panel_cmp,
  post_cmp,
  by = c(
    "N",
    "rho_y",
    "treat_prob",
    "flood_exposure_rate",
    "tau_flood",
    "tau_interaction"
  )
) |>
  mutate(gain_panel_minus_post = panel_power - post_power)

summary_tbl <- tibble(
  metric = c(
    "flood_post_median_power",
    "flood_panel_median_power",
    "flood_panel_minus_post_median_gain",
    "flood_panel_minus_post_gain_q10",
    "flood_panel_minus_post_gain_q90",
    "flood_post_share_ge_80_power",
    "flood_panel_share_ge_80_power",
    "interaction_post_median_power",
    "interaction_panel_median_power",
    "interaction_post_share_ge_80_power",
    "interaction_panel_share_ge_80_power"
  ),
  value = c(
    median(post_cmp$post_power),
    median(panel_cmp$panel_power),
    median(joined$gain_panel_minus_post),
    as.numeric(quantile(joined$gain_panel_minus_post, 0.10)),
    as.numeric(quantile(joined$gain_panel_minus_post, 0.90)),
    mean(post_cmp$post_power >= 0.80),
    mean(panel_cmp$panel_power >= 0.80),
    median(inter_power$power[!inter_power$panel]),
    median(inter_power$power[inter_power$panel]),
    mean(inter_power$power[!inter_power$panel] >= 0.80),
    mean(inter_power$power[inter_power$panel] >= 0.80)
  )
)

by_attrition_gain <- joined |>
  group_by(attrition_rate) |>
  summarise(median_gain = median(gain_panel_minus_post), .groups = "drop")

by_rho_gain <- joined |>
  group_by(rho_y) |>
  summarise(median_gain = median(gain_panel_minus_post), .groups = "drop")

by_floodrate_panel_floodpower <- panel_cmp |>
  group_by(flood_exposure_rate) |>
  summarise(median_power = median(panel_power), .groups = "drop")

by_treatprob_panel_interpower <- inter_power |>
  filter(panel) |>
  group_by(treat_prob) |>
  summarise(median_power = median(power), .groups = "drop")

snapshot <- list(
  run_meta = list(
    timestamp = Sys.time(),
    sims = sims,
    grids = list(
      N = Ns,
      rho_y = rhos,
      attrition_rate_panel = attrs,
      treat_prob = treats,
      flood_exposure_rate = floods,
      tau_flood = tau_floods,
      tau_interaction = tau_inters
    )
  ),
  summary = summary_tbl,
  tables = list(
    by_attrition_gain = by_attrition_gain,
    by_rho_gain = by_rho_gain,
    by_floodrate_panel_floodpower = by_floodrate_panel_floodpower,
    by_treatprob_panel_interpower = by_treatprob_panel_interpower
  )
)

# -- write outputs -------------------------------------------------------------
saveRDS(snapshot, file.path(out_dir, "quick_findings_snapshot.rds"))
write.csv(summary_tbl, file.path(out_dir, "quick_findings_summary.csv"), row.names = FALSE)
write.csv(by_attrition_gain, file.path(out_dir, "quick_findings_by_attrition_gain.csv"), row.names = FALSE)
write.csv(by_rho_gain, file.path(out_dir, "quick_findings_by_rho_gain.csv"), row.names = FALSE)
write.csv(by_floodrate_panel_floodpower, file.path(out_dir, "quick_findings_by_floodrate_panel_floodpower.csv"), row.names = FALSE)
write.csv(by_treatprob_panel_interpower, file.path(out_dir, "quick_findings_by_treatprob_panel_interpower.csv"), row.names = FALSE)

# -- console summary -----------------------------------------------------------
message("\nQuick findings snapshot written to: ", out_dir)
message("- quick_findings_snapshot.rds")
message("- quick_findings_summary.csv")
message("- quick_findings_by_attrition_gain.csv")
message("- quick_findings_by_rho_gain.csv")
message("- quick_findings_by_floodrate_panel_floodpower.csv")
message("- quick_findings_by_treatprob_panel_interpower.csv")

message("\nHeadline metrics:")
print(summary_tbl)
