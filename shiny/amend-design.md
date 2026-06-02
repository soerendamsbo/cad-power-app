# How to amend the WP2 power analysis design

All design logic lives in one file: `shiny/design_builder.R`.
After editing it, re-run `Rscript shiny/precompute.R` from the
project root, then relaunch the app.

---

## Adding a new estimand (effect)

An estimand is a target quantity — what the study is trying to estimate.

**Step 1 — Add potential outcomes to the model** (`declare_model` block,
~line 82). Define the new outcome's potential-outcome columns. Example —
adding a democratic-trust outcome `Y_trust`:

```r
Y_trust_Z_0 = rho_y * Y_pre_latent + sqrt(1 - rho_y^2) * rnorm(N),
Y_trust_Z_1 = Y_trust_Z_0 + tau,   # same tau for now; give it own parameter if needed
```

**Step 2 — Reveal the outcome in measurement** (`declare_measurement`
block, ~line 159):

```r
Y_trust = ifelse(Z == 1, Y_trust_Z_1, Y_trust_Z_0)
```

**Step 3 — Declare the estimand** (`declare_inquiry` block, ~line 149,
marked `## ADD NEW ESTIMANDS HERE`):

```r
ATE_trust = mean(Y_trust_Z_1 - Y_trust_Z_0)
```

**Step 4 — Add estimators** (estimator blocks, ~line 182, marked
`## ADD NEW ESTIMATORS HERE`). One per design type that should cover the
new estimand:

```r
declare_estimator(Y_trust ~ Z,
  .method = lm, inquiry = "ATE_trust", term = "Z",
  label = "post_only_trust")
```

Add the corresponding panel ANCOVA version if needed:

```r
declare_estimator(Y_trust ~ Z + Y_pre_latent,
  .method = lm, inquiry = "ATE_trust", term = "Z",
  label = "panel_adjusted_trust")
```

**Step 5 — Register the display label** (`ESTIMATOR_LABELS` vector,
~line 201, marked `## ADD NEW ESTIMATOR LABELS HERE`):

```r
post_only_trust       = "Post-only · Trust",
panel_adjusted_trust  = "Panel adjusted · Trust"
```

That's all. The app's plots, tables, and MDE tab pick up new estimators
automatically from the data.

---

## Adding a new design type

A design type is an assembly of model + inquiry + data strategy +
estimator (e.g. a heteroskedasticity-robust variant, a subgroup analysis,
or a two-stage design).

1. Add a new `else if (type == "my_type")` branch in the assembly block
   at the bottom of `build_wp2_design()` (~line 219, marked
   `## ADD NEW DESIGN TYPES HERE`).
2. Register it in `DESIGN_REGISTRY` (~line 186):
   ```r
   list(type = "my_type", label = "My design label",
        uses_attrition = TRUE)   # FALSE if attrition params don't apply
   ```
3. Add estimator labels in `ESTIMATOR_LABELS` as above.

---

## Changing the DGP assumptions

All population parameters (demographic proportions, ideology–climate
correlation, baseline outcome loadings, heterogeneous-effect coefficients)
are in the `declare_model()` call (~line 82). Edit them directly.

The cross-wave stability parameter `rho_y` is a function argument and is
already part of the precompute grid — no further changes needed if you
only want to vary it.

---

## Expanding the parameter grid

Edit the grid vectors at the top of `precompute.R` (clearly marked with
`## ADD:` and `## RAISE:` comments). Each extra value multiplies total
cells and runtime proportionally; see the per-value cost annotations in
that file.

---

## Files not to touch for design changes

- `app.R` — only needs changes for UI layout or new filter controls.
  New estimators appear in existing plots automatically.
- `precompute.R` — only needs changes for grid coverage.
  Design logic must stay in `design_builder.R`.
