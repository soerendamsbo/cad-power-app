# Findings Grid: Run Plan and Output Instructions

## Step 1 — Run the simulation grid

From the project root:

```
Rscript shiny/findings_grid.R
```

This runs ~4,806 design cells (WP3 post-only + WP3 panel) × 200 simulations each
using `shiny/design_builder.R`. Parallel execution via `{future}` if available.
Expected runtime: 30–60 minutes on a modern laptop with parallel.

Output written to `shiny/findings_grid_out/`:
- `findings_grid.csv`  — one row per (cell × estimator × inquiry), ~77K rows
- `findings_grid.rds`  — list(results = <data.frame>, meta = <list>)

## Step 2 — Analyse the output

Load with:

```r
library(dplyr)
r <- readRDS(here::here("shiny", "findings_grid_out", "findings_grid.rds"))$results
```

Key columns: `design_type`, `N`, `tau`, `rho_y`, `attrition_rate`,
`flood_exposure_rate`, `tau_flood`, `tau_interaction`, `flood_response_boost`,
`treat_prob`, `estimator`, `inquiry`, `term`, `outcome_scale`,
`power`, `bias`, `rmse`, `coverage`, `mean_se`, `type_s_error`,
`n_sims`, `expected_flooded_n`.

Focus on `outcome_scale == "Latent SD"` rows (filter `!grepl("1_5", estimator)`).

Key inquiry values:
- `flood_avg_latent`              — flood main effect (primary bottleneck)
- `flood_survey_contrast_latent`  — interaction effect (hardest to power)
- `survey_avg_latent`             — survey experiment main effect

## Step 3 — Derive findings to update app.R

Compute and summarise these specific questions:

**Q1. N × flood_exposure_rate frontier for flood-effect power**
Filter: `inquiry == "flood_avg_latent"`, `design_type == "wp3_post_only"`, `tau_flood == -0.10`
Group by N, flood_exposure_rate. Report where median power crosses 0.80.

**Q2. Panel vs post-only gain**
Filter: `inquiry == "flood_avg_latent"`, latent scale.
Join panel and post-only on shared parameters (N, tau, flood_exposure_rate, tau_flood,
tau_interaction). Compute `gain = panel_power - post_power`.
Summarise median gain, q10/q90, and cross-tab by rho_y and attrition_rate.

**Q3. Interaction power**
Filter: `inquiry == "flood_survey_contrast_latent"`, latent scale.
Report median power by design_type, and share of cells reaching 0.80.
Does interaction power ever reach 0.80 in the grid?

**Q4. flood_response_boost sensitivity**
Filter panel rows: `inquiry == "flood_avg_latent"`, compare `flood_response_boost == 0`
vs `flood_response_boost == 0.05`. Summarise whether retention bias
meaningfully shifts power or bias.

**Q5. Survey-effect power (tau sensitivity)**
Filter: `inquiry == "survey_avg_latent"`, `design_type == "wp3_post_only"`.
Report power by tau (0.05, 0.10, 0.20) and N.

## Step 4 — Update app.R

The section to update is `preliminary_findings_content` in `shiny/app.R`,
starting around line 1021. It is a `tagList(...)` with hardcoded bullet points.

Replace the bullet-point text under each `h4()` heading with findings derived
from Step 3. Keep the same HTML structure (tagList, h4, tags$ul, tags$li,
tags$p). Update the caveat paragraph (last `tags$p`) to reflect the new grid
parameters and sims=200.

Do NOT change the structure of the app, add new tabs, or touch any other
section of app.R.
