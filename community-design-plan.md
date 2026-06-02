# Implementation Plan: WP3 Community + Individual Flooding Design

## Rationale and key decisions

**Why a separate builder function, not extending `build_wp2_design`.**
The community design has a structurally different parameter set (J, p_comm_flooded, p_indiv_given_comm vs. flood_exposure_rate), a different DGP (hierarchical, not IID), and different estimands. Extending the existing function's `type` argument would require conditionally ignoring half the parameters in each branch, which obscures both designs.

**Why community-level cluster-robust SEs across the board, not community FE.**
Community FE would absorb `community_flooded` entirely — it is constant within each community by construction, so FE collinears it out and the community effect cannot be estimated. We use `lm_robust(..., clusters = community_id, se_type = "CR2")` without community FE. This estimates both effects in one model with conservative, valid inference for both.

**Why the additive formula over a categorical `exposure_status` encoding.**
`Y ~ Z * community_flooded + Z * individual_flooded` maps one-to-one to parameters (tau_community, tau_individual, tau), avoids dummy-variable encoding choices, and makes the structural-zero constraint (can't have individual_flooded=1 and community_flooded=0) implicit in the DGP rather than requiring re-parameterisation.

**Why no survey × flood interaction in this design.**
The existing WP3 design has `tau_interaction` (Z × flooded). Keeping it in the community design means tracking two additional interaction terms (Z×community, Z×individual), which quickly overwhelms the power budget and obscures the main diagnostic. Drop interactions in the community DGP for this version; the Methods section should note this assumption.

**Why latent scale only (no Likert) for the community design.**
Likert adds parallel potential outcomes, measurement blocks, and estimator branches. For a power analysis focused on the effective-N collapse story, the latent scale is sufficient. Add Likert later if needed.

**The critical power message the app must surface.**
The community effect is identified between communities; effective N ≈ J × p_comm_flooded (number of flooded communities). The individual/property effect is identified within flooded communities; effective N ≈ N × p_comm_flooded × p_indiv_given_comm (number of directly flooded people). These can differ by an order of magnitude. The diagnosand table will surface this automatically via `mean_se`.

---

## Files to change

1. `shiny/design_builder.R`
2. `shiny/precompute.R`
3. `shiny/app.R`

## Files that do not change

All existing WP2 and WP3 designs and their precomputed results are preserved unchanged. New rows are appended to `results.rds`. New columns in those rows (J, p_comm_flooded, etc.) are NA for existing rows.

---

## Step 1 — `design_builder.R`

### 1a. New estimator function

Add `wp3_community_effects_estimator()` after the existing `wp3_effects_estimator()`. Signature:

```r
wp3_community_effects_estimator <- function(
  data,
  outcome,
  adjust_pre = FALSE,
  term = NULL
)
```

- Fits `estimatr::lm_robust(Y ~ Z * community_flooded + Z * individual_flooded [+ Y_pre_latent], data, clusters = community_id, se_type = "CR2")`.
- Returns a data frame with columns `term, estimate, std.error, statistic, p.value, conf.low, conf.high` for the terms: `"Z"`, `"community_flooded"`, `"individual_flooded"`, `"total_direct_effect_latent"`, `"Z:community_flooded"`, `"Z:individual_flooded"`.
- Extract coefficient terms directly from the `lm_robust` fit coefficients/SEs/CIs using `tidy(fit)` or manual extraction.
- Compute `total_direct_effect_latent` as the linear combination `community_flooded + individual_flooded`, using the cluster-robust variance-covariance matrix from the same `lm_robust` fit.

### 1b. New builder function

Add `build_wp3_community_design()` with parameters:

| Parameter | Default | Description |
|---|---|---|
| `N` | — | Total sample size |
| `J` | 200 | Number of communities |
| `p_comm_flooded` | 0.15 | Proportion of communities that flood |
| `p_indiv_given_comm` | 0.40 | Proportion with direct property damage within flooded communities |
| `tau` | 0.10 | Survey experiment effect |
| `tau_community` | -0.10 | Community-level flooding effect |
| `tau_individual` | -0.05 | Incremental individual/property effect within flooded communities |
| `rho_y` | 0.40 | Cross-wave correlation |
| `attrition_rate` | 0.30 | Baseline wave-2 dropout |
| `differential_attrition` | 0.00 | Log-odds dropout penalty for survey-treated |
| `treat_prob` | 0.50 | Survey experiment assignment probability |
| `flood_response_boost` | 0.00 | Log-odds retention boost for community-flooded |
| `type` | — | `"wp3_community_post_only"` or `"wp3_community_panel"` |

**DGP (`declare_model`):**

```r
community_id = sample(seq_len(J), N, replace = TRUE),
community_flooded_j = rbinom(J, 1, p_comm_flooded),
community_flooded = community_flooded_j[community_id],
individual_flooded = community_flooded * rbinom(N, 1, p_indiv_given_comm),
# ... standard demographics, Y_pre_latent, Y0_latent, tau_z_i as in existing builder ...
# Potential outcomes (6-cell, no Z x flood interaction):
Y_latent_Z0_C0    = Y0_latent,
Y_latent_Z1_C0    = Y0_latent + tau_z_i,
Y_latent_Z0_C1_I0 = Y0_latent + tau_community,
Y_latent_Z1_C1_I0 = Y0_latent + tau_z_i + tau_community,
Y_latent_Z0_C1_I1 = Y0_latent + tau_community + tau_individual,
Y_latent_Z1_C1_I1 = Y0_latent + tau_z_i + tau_community + tau_individual
```

**Inquiries (`declare_inquiry`):**

```r
community_effect_latent    = mean(Y_latent_Z0_C1_I0 - Y_latent_Z0_C0),
individual_effect_latent   = mean(Y_latent_Z0_C1_I1 - Y_latent_Z0_C1_I0),
survey_effect_latent       = mean(Y_latent_Z1_C0 - Y_latent_Z0_C0),
total_direct_effect_latent = mean(Y_latent_Z0_C1_I1 - Y_latent_Z0_C0)
```

Note: `individual_effect_latent` is the incremental direct property-damage effect among people in flooded communities. Under the additive DGP above, it equals `tau_individual`, not a prevalence-weighted population burden. This is consistent with how `individual_flooded` enters the regression coefficient. Note this in the Methods section.

**Measurement:** Reveal observed `Y_latent` from the 6-cell potential outcome lookup using `community_flooded` and `individual_flooded`. No Likert variant for now.

**Attrition:** Same logistic model as existing, with `community_flooded` in place of `flooded`. Use parameter name `flood_response_boost` for continuity.

**Estimators:** Two estimator declarations per design type (post-only and panel):

```r
declare_estimator(
  handler = label_estimator(wp3_community_effects_estimator),
  outcome = "Y_latent",
  adjust_pre = FALSE,   # TRUE for panel
  inquiry = c(
    "survey_effect_latent",
    "community_effect_latent",
    "individual_effect_latent",
    "total_direct_effect_latent",
    NA_character_,
    NA_character_
  ),
  term = c(
    "Z",
    "community_flooded",
    "individual_flooded",
    "total_direct_effect_latent",
    "Z:community_flooded",
    "Z:individual_flooded"
  ),
  label = "wp3_community_post_latent"  # or "_panel_"
)
```

Map `term` to `inquiry` so DeclareDesign correctly links each returned row to its estimand. `term` and `inquiry` must have the same length when passed through `label_estimator()`. The terms `Z:community_flooded` and `Z:individual_flooded` do not have matching inquiries in this version and should be assigned `NA_character_`.

**Assembly:** Two branches (`wp3_community_post_only`, `wp3_community_panel`) analogous to existing WP3 branches.

### 1c. Registry and labels

Append to `DESIGN_REGISTRY`:

```r
list(type = "wp3_community_post_only",
     label = "WP3 - Community + property flooding (post-only)",
     uses_attrition = FALSE),
list(type = "wp3_community_panel",
     label = "WP3 - Community + property flooding (panel)",
     uses_attrition = TRUE)
```

Append to `ESTIMATOR_LABELS`:

```r
wp3_community_post_latent  = "WP3 Community - Post-only (latent SD)",
wp3_community_panel_latent = "WP3 Community - Panel (latent SD)"
```

Add two new exported constants at the bottom of the file:

```r
WP3_COMMUNITY_EFFECT_ORDER <- c(
  "community_flooded",
  "individual_flooded",
  "Z",
  "total_direct_effect_latent",
  "Z:community_flooded",
  "Z:individual_flooded"
)

WP3_COMMUNITY_EFFECT_LABELS <- c(
  community_flooded             = "Community flooding effect",
  individual_flooded            = "Individual/property effect (within flooded communities)",
  Z                             = "Survey experiment effect",
  total_direct_effect_latent    = "Total direct-damage effect",
  `Z:community_flooded`         = "Survey x community interaction",
  `Z:individual_flooded`        = "Survey x individual interaction"
)

WP3_COMMUNITY_EFFECT_X <- c(
  community_flooded          = "tau_community",
  individual_flooded         = "tau_individual",
  Z                          = "tau",
  total_direct_effect_latent = "tau_individual",
  `Z:community_flooded`      = "tau",
  `Z:individual_flooded`     = "tau"
)

WP3_COMMUNITY_EFFECT_X_LABELS <- c(
  tau = "Survey experiment effect tau (SD units)",
  tau_community = "Community flooding effect tau_community (SD units)",
  tau_individual = "Individual/property effect tau_individual (SD units)"
)
```

---

## Step 2 — `precompute.R`

### 2a. New parameter grids

Add at the top of the parameter grids section:

```r
# WP3 community design grids
J_grid                <- c(100, 200)
p_comm_flooded_grid   <- c(0.10, 0.20)
p_indiv_grid          <- c(0.40)
tau_community_grid    <- c(-0.10, -0.05)
tau_individual_grid   <- c(-0.05)
# tau uses tau_grid_wp3 (already defined)
```

Keep this initial precompute grid deliberately lean because the community estimator uses CR2 cluster-robust inference, which is materially slower than the existing IID OLS estimators. Expand these vectors only for a final precompute or use the custom simulation tab for targeted scenarios.

### 2b. New design grids

```r
grid_wp3_community_post <- expand.grid(
  N                  = N_grid,
  J                  = J_grid,
  tau                = tau_grid_wp3,
  rho_y              = rho_grid,
  p_comm_flooded     = p_comm_flooded_grid,
  p_indiv_given_comm = p_indiv_grid,
  tau_community      = tau_community_grid,
  tau_individual     = tau_individual_grid,
  treat_prob         = treat_prob_grid,
  stringsAsFactors   = FALSE
)

grid_wp3_community_panel <- expand.grid(
  N                      = N_grid,
  J                      = J_grid,
  tau                    = tau_grid_wp3,
  rho_y                  = rho_grid,
  p_comm_flooded         = p_comm_flooded_grid,
  p_indiv_given_comm     = p_indiv_grid,
  tau_community          = tau_community_grid,
  tau_individual         = tau_individual_grid,
  treat_prob             = treat_prob_grid,
  attrition_rate         = wp3_attrition_grid,
  differential_attrition = wp3_differential_attrition_grid,
  flood_response_boost   = flood_response_boost_grid,
  stringsAsFactors       = FALSE
)
```

### 2c. Build and diagnose

Add build/diagnose loops for `wp3_community_post_designs` and `wp3_community_panel_designs` using `build_wp3_community_design()`, following the same pattern as the existing WP3 loops. Pass `tau = grid_wp3_community_post$tau[i]` or `tau = grid_wp3_community_panel$tau[i]` to the builder and keep the joined results column named `tau` so existing app paths continue to work.

### 2d. Join parameters and bind

Join parameters back to diagnosands using the same `parse_idx` pattern. NA-fill the community-specific columns for non-community rows and vice versa. Add `results_wp3_community_post` and `results_wp3_community_panel` to the `bind_rows` call.

Update `attr(results, "grid_meta")` to include the new grid vectors:

```r
J_grid                = J_grid,
p_comm_flooded_grid   = p_comm_flooded_grid,
p_indiv_grid          = p_indiv_grid,
tau_community_grid    = tau_community_grid,
tau_individual_grid   = tau_individual_grid,
```

---

## Step 3 — `app.R`

### 3a. Grid metadata extraction

Add the community columns to the backward-compatibility block that currently creates missing legacy columns after `results <- readRDS("results.rds")`:

```r
for (nm in c(
  "J",
  "p_comm_flooded",
  "p_indiv_given_comm",
  "tau_community",
  "tau_individual"
)) {
  if (!nm %in% names(results)) {
    results[[nm]] <- NA_real_
  }
}
```

After the existing `meta_or_col` calls, add:

```r
J_vals              <- meta_or_col("J_grid", "J", default = c(100, 200))
p_comm_flooded_vals <- meta_or_col("p_comm_flooded_grid", "p_comm_flooded", default = c(0.10, 0.20))
p_indiv_vals        <- meta_or_col("p_indiv_grid", "p_indiv_given_comm", default = c(0.40))
tau_community_vals  <- meta_or_col("tau_community_grid", "tau_community", default = c(-0.10, -0.05))
tau_individual_vals <- meta_or_col("tau_individual_grid", "tau_individual", default = c(-0.05))
```

### 3b. Precomputed sidebar

Add five new `selectInput` / `conditionalPanel` controls for the community parameters, following the same style as existing WP3 controls (`filter_J`, `filter_p_comm_flooded`, `filter_p_indiv`, `filter_tau_community`, `filter_tau_individual`). All five should be conditionally hidden when faceting by the same variable.

Add the five new parameters to `FACET_CHOICES` and `FACET_LABELS`.

Update `filter_subtitle()` to include the new parameters.

### 3c. New community design plot section in dashboard

Add a new `nav_panel("WP3 Community Design")` in `dashboard_main`. Inside it:

- A note explaining: *"The community effect is identified between communities (effective N ≈ flooded communities). The individual/property effect is identified within flooded communities. Mean SE in the numerical table shows this difference concretely."*
- A plot grid for the main estimands: community effect, individual/property effect, survey effect, total direct-damage effect.

Add a `wp3_community_effect_grid_ui()` helper analogous to `wp3_effect_grid_ui()` but using `WP3_COMMUNITY_EFFECT_ORDER`.

Register community effect plots via a new `register_wp3_community_effect_plots()` call, using a `wp3_community_effect_plot()` function that:
- Filters `results` to `design_type %in% c("wp3_community_post_only", "wp3_community_panel")`
- Applies the community-specific sidebar filters (J, p_comm_flooded, p_indiv, tau_community, tau_individual) in addition to the shared filters (N, rho_y, treat_prob)
- Uses `WP3_COMMUNITY_EFFECT_X` to choose the x-axis for each effect. Hold `tau_community` fixed when plotting against `tau_individual`, hold `tau_individual` fixed when plotting against `tau_community`, and hold both fixed for survey-effect plots.

### 3d. Custom simulation tab

Add `wp3_community_post_only` and `wp3_community_panel` to `CUSTOM_WAVE_CHOICES` (or a separate checkbox group "Community design modes").

Add `numericInput` controls for the five community parameters inside the custom simulation card, wrapped in a `conditionalPanel` that shows only when a community design type is selected:

- `custom_J` (Number of communities, default 200)
- `custom_p_comm_flooded` (Share of flooded communities, default 0.15)
- `custom_p_indiv` (Direct damage rate within flooded communities, default 0.40)
- `custom_tau_community` (Community flooding effect, default -0.10)
- `custom_tau_individual` (Individual/property effect, default -0.05)

In the `eventReactive(input$run_custom)` block: for community design types, call `build_wp3_community_design()` instead of `build_wp2_design()`, passing the community-specific inputs.

Add a conditional plot section in the custom results that shows community estimand plots when community design results are present, separate from the existing WP3 effect grid.

Update `custom_memory_columns()` to include the new parameters. Update `custom_memory_rows()` to record community-design power estimates. New columns are NA for non-community runs.

### 3e. Filter reactive and cell count

Update the `plot_data` reactive to pass through the new community parameter filters (with `is.na()` passthrough for non-community rows, same as existing WP3 parameters).

Update `n_cells` calculation to include the new community design columns.

Update the numerical results and custom-results tables to include `J`, `p_comm_flooded`, `p_indiv_given_comm`, `tau_community`, and `tau_individual`, with `NA` values displayed for non-community designs.

### 3f. Methods section

Add a paragraph to `methods_content` under "Model" describing the community design DGP:

- Hierarchical assignment: N individuals assigned to J communities, community flood status drawn as J Bernoullis with probability `p_comm_flooded`, individual property exposure drawn conditionally on community flood status.
- Additive potential outcomes without nonzero Z × flood effects; the regression still includes `Z:community_flooded` and `Z:individual_flooded` diagnostic terms, but they have no nonzero target estimand in this version.
- Estimator: `lm_robust(Y ~ Z * community_flooded + Z * individual_flooded [+ Y_pre_latent], clusters = community_id, se_type = "CR2")`.
- IID standard errors would be invalid here; cluster-robust SEs at the community level are used for both effects.
- Effective N for the community effect depends on the number of flooded communities, not the individual N.
- The individual/property coefficient targets the incremental direct-damage effect conditional on community flooding (`tau_individual`), not a prevalence-weighted population-average burden.

---

## Verification checkpoints

**After Step 1:** `build_wp3_community_design(N=500, J=100, p_comm_flooded=0.20, p_indiv_given_comm=0.40, tau=0.10, tau_community=-0.10, tau_individual=-0.05, rho_y=0.40, type="wp3_community_post_only")` runs without error and `diagnose_design(...)` returns diagnosands for all four inquiries with plausible values.

**After Step 2:** `Rscript precompute.R` completes and `results.rds` contains rows for `design_type %in% c("wp3_community_post_only", "wp3_community_panel")` with non-NA values in `J`, `p_comm_flooded`, `p_indiv_given_comm`, `tau_community`, `tau_individual`.

**After Step 3:** App launches without error; the "WP3 Community Design" tab shows plots; the custom simulation runs a community design and produces power estimates for `community_flooded` and `individual_flooded` terms; the sidebar community filters filter results correctly for the precomputed dashboard.
