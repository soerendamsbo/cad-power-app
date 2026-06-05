# app.R - WP2/WP3 Power Analysis Explorer
#
# Launch: shiny::runApp("app.R")
# Requires results.rds (produced by precompute.R) in the same directory.
#
# Required packages:
#   install.packages(c("shiny", "bslib", "dplyr", "ggplot2", "DT", "scales"))
#
# ── Sections ──────────────────────────────────────────────────────────────────
#   1. Libraries and setup
#   2. Load pre-computed data
#   3. Static text (methods, interpretation)
#   4. UI
#   5. Server
#   6. shinyApp()
# ─────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# 1. Libraries and setup
# ══════════════════════════════════════════════════════════════════════════════

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(DT)
library(scales)

source("design_builder.R") # latent_to_likert, diagnosands_power,
# build_wp2_design, DESIGN_REGISTRY, ESTIMATOR_LABELS

# ══════════════════════════════════════════════════════════════════════════════
# 2. Load pre-computed data
# ══════════════════════════════════════════════════════════════════════════════

if (!file.exists("results.rds")) {
  stop(
    "results.rds not found.\n",
    "Run precompute.R first:\n",
    "  Rscript shiny/precompute.R\n",
    "Then relaunch the app."
  )
}

results <- readRDS("results.rds")

# Backward compatibility for older precompute outputs.
results <- results |>
  mutate(
    design_type = dplyr::recode(
      design_type,
      post_only = "wp2_post_only",
      panel = "wp2_panel"
    ),
    outcome_scale = gsub("Likert 1.?5", "Likert 1-5", outcome_scale)
  )

for (nm in c(
  "flood_exposure_rate",
  "treat_prob",
  "tau_flood",
  "tau_interaction",
  "flood_response_boost",
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

grid_meta <- attr(results, "grid_meta")

col_vals <- function(col, default = numeric(0)) {
  if (!col %in% names(results)) {
    return(default)
  }
  vals <- sort(unique(stats::na.omit(results[[col]])))
  if (length(vals) == 0) default else vals
}

meta_or_col <- function(meta_name, col_name, default = numeric(0)) {
  vals <- grid_meta[[meta_name]]
  if (is.null(vals)) {
    return(col_vals(col_name, default))
  }
  vals <- sort(vals)
  if (length(vals) == 0) col_vals(col_name, default) else vals
}

# Extract grid values for UI controls (from metadata, not re-derived from data)
N_vals <- meta_or_col("N_grid", "N")
tau_vals <- meta_or_col("tau_grid", "tau")
rho_vals <- meta_or_col("rho_grid", "rho_y")
attr_vals <- meta_or_col("attrition_grid", "attrition_rate", default = c(0.30))
diff_vals <- meta_or_col(
  "differential_attrition_grid",
  "differential_attrition",
  default = c(0.00)
)
flood_vals <- meta_or_col(
  "flood_exposure_grid",
  "flood_exposure_rate",
  default = c(0.10, 0.20)
)
treat_vals <- meta_or_col(
  "treat_prob_grid",
  "treat_prob",
  default = c(0.50, 0.65)
)
tau_flood_vals <- meta_or_col(
  "tau_flood_grid",
  "tau_flood",
  default = c(-0.10, -0.05)
)
tau_int_vals <- meta_or_col(
  "tau_interaction_grid",
  "tau_interaction",
  default = c(-0.08, -0.04, 0.04, 0.08)
)
flood_boost_vals <- meta_or_col(
  "flood_response_boost_grid",
  "flood_response_boost",
  default = c(0.05, 0.10)
)
J_vals <- meta_or_col("J_grid", "J", default = c(100, 200))
p_comm_flooded_vals <- meta_or_col(
  "p_comm_flooded_grid",
  "p_comm_flooded",
  default = c(0.10, 0.20)
)
p_indiv_vals <- meta_or_col(
  "p_indiv_grid",
  "p_indiv_given_comm",
  default = c(0.40)
)
tau_community_vals <- meta_or_col(
  "tau_community_grid",
  "tau_community",
  default = c(-0.10, -0.05)
)
tau_individual_vals <- meta_or_col(
  "tau_individual_grid",
  "tau_individual",
  default = c(-0.05)
)
sims_used <- if (!is.null(grid_meta$sims)) grid_meta$sims else NA_integer_
precompute_ts <- if (!is.null(grid_meta$timestamp)) {
  format(grid_meta$timestamp, "%Y-%m-%d %H:%M")
} else {
  "unknown"
}
n_cells <- results |>
  distinct(
    design_type,
    N,
    tau,
    rho_y,
    attrition_rate,
    differential_attrition,
    flood_exposure_rate,
    treat_prob,
    tau_flood,
    tau_interaction,
    flood_response_boost,
    J,
    p_comm_flooded,
    p_indiv_given_comm,
    tau_community,
    tau_individual
  ) |>
  nrow()

DESIGN_TYPES <- vapply(DESIGN_REGISTRY, function(x) x$type, character(1))
DESIGN_LABELS <- vapply(DESIGN_REGISTRY, function(x) x$label, character(1))
DESIGN_CHOICES <- setNames(DESIGN_TYPES, DESIGN_LABELS)
DESIGN_LABEL_LOOKUP <- setNames(DESIGN_LABELS, DESIGN_TYPES)
CUSTOM_WAVE_CHOICES <- c(
  "Waves = 1 (post only)" = "post_only",
  "Waves = 2 (pre-post)" = "panel"
)
CUSTOM_DESIGN_TYPE_ORDER <- c(
  "wp3_community_post_only",
  "wp3_community_panel"
)

# ── Plotting helpers ──────────────────────────────────────────────────────────

FACET_CHOICES <- c(
  "Sample size (N)" = "N",
  "Cross-wave correlation (ρ)" = "rho_y",
  "Attrition rate" = "attrition_rate",
  "Survey experiment retention penalty" = "differential_attrition",
  "Flood exposure rate" = "flood_exposure_rate",
  "Survey experiment assignment probability" = "treat_prob",
  "Flooding effect" = "tau_flood",
  "Flood x Survey interaction" = "tau_interaction",
  "Flood retention boost" = "flood_response_boost",
  "Number of communities (J)" = "J",
  "Flooded community share" = "p_comm_flooded",
  "Individual flooding rate within flooded communities" = "p_indiv_given_comm",
  "Community flooding effect" = "tau_community",
  "Individual flooding effect" = "tau_individual"
)

FACET_LABELS <- c(
  N = "N",
  rho_y = "\u03c1 (cross-wave correlation)",
  attrition_rate = "Attrition rate",
  differential_attrition = "Survey experiment retention penalty",
  flood_exposure_rate = "Flood exposure rate",
  treat_prob = "Survey experiment assignment probability",
  tau_flood = "Flooding effect",
  tau_interaction = "Flood x Survey interaction",
  flood_response_boost = "Flood retention boost",
  J = "J",
  p_comm_flooded = "Flooded community share",
  p_indiv_given_comm = "Individual flooding rate",
  tau_community = "Community flooding effect",
  tau_individual = "Individual flooding effect"
)

WP3_EFFECT_LABELS <- c(
  survey_avg = "Average survey experiment effect",
  survey_when_not_flooded = "Survey experiment effect: not flooded",
  survey_when_flooded = "Survey experiment effect: flooded",
  survey_flood_contrast = "Survey-effect contrast",
  flood_avg = "Average flooding effect",
  flood_when_survey_control = "Flooding effect: survey control",
  flood_when_survey_treated = "Flooding effect: survey treated",
  flood_survey_contrast = "Flooding-effect contrast"
)

WP3_EFFECT_X <- c(
  survey_avg = "tau",
  survey_when_not_flooded = "tau",
  survey_when_flooded = "tau",
  survey_flood_contrast = "tau_interaction",
  flood_avg = "tau_flood",
  flood_when_survey_control = "tau_flood",
  flood_when_survey_treated = "tau_flood",
  flood_survey_contrast = "tau_interaction"
)

WP3_EFFECT_X_LABELS <- c(
  tau = "Survey experiment effect tau (SD units)",
  tau_flood = "Flood exposure effect tau_flood (SD units)",
  tau_interaction = "Survey x flood interaction (SD units)"
)

WP3_EFFECT_ORDER <- c(
  "survey_avg",
  "flood_avg",
  "survey_when_not_flooded",
  "survey_when_flooded",
  "flood_when_survey_control",
  "flood_when_survey_treated",
  "survey_flood_contrast",
  "flood_survey_contrast"
)

wp3_effect_grid_ui <- function(prefix, height = "300px") {
  tagList(
    h4("Average effects"),
    fluidRow(
      column(6, plotOutput(paste0(prefix, "_survey_avg"), height = height)),
      column(6, plotOutput(paste0(prefix, "_flood_avg"), height = height))
    ),
    tags$hr(),
    h4("Conditional effects"),
    h5("Survey experiment effects"),
    fluidRow(
      column(
        6,
        plotOutput(paste0(prefix, "_survey_when_not_flooded"), height = height)
      ),
      column(
        6,
        plotOutput(paste0(prefix, "_survey_when_flooded"), height = height)
      )
    ),
    h5("Flooding effects"),
    fluidRow(
      column(
        6,
        plotOutput(
          paste0(prefix, "_flood_when_survey_control"),
          height = height
        )
      ),
      column(
        6,
        plotOutput(
          paste0(prefix, "_flood_when_survey_treated"),
          height = height
        )
      )
    ),
    tags$hr(),
    h4("Contrasts"),
    fluidRow(
      column(
        6,
        plotOutput(paste0(prefix, "_survey_flood_contrast"), height = height)
      ),
      column(
        6,
        plotOutput(paste0(prefix, "_flood_survey_contrast"), height = height)
      )
    )
  )
}

effect_output_id <- function(prefix, effect_term) {
  paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", effect_term))
}

wp3_community_effect_grid_ui <- function(prefix, height = "300px") {
  tagList(
    h4("Primary effects"),
    fluidRow(
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "Z"),
          height = height
        )
      ),
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "community_flooded"),
          height = height
        )
      ),
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "individual_flooded"),
          height = height
        )
      )
    ),
    tags$hr(),
    h4("Combined and diagnostic effects"),
    fluidRow(
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "total_direct_effect_latent"),
          height = height
        )
      ),
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "Z:community_flooded"),
          height = height
        )
      ),
      column(
        4,
        plotOutput(
          effect_output_id(prefix, "Z:individual_flooded"),
          height = height
        )
      )
    )
  )
}

# ── Plot helpers ──────────────────────────────────────────────────────────────

# Wrap long estimator labels for plot legends.
wrap_label <- function(x, width = 28) {
  vapply(
    x,
    function(s) paste(strwrap(s, width), collapse = "\n"),
    character(1),
    USE.NAMES = FALSE
  )
}

param_note <- function(...) {
  div(class = "form-text text-muted mb-2", ...)
}

is_valid_number <- function(x) {
  length(x) == 1 && is.numeric(x) && is.finite(x)
}

effect_size_helper_value <- function(
  raw_effect,
  outcome_sd,
  scale_min,
  scale_max
) {
  if (!is_valid_number(raw_effect)) {
    return(NA_real_)
  }

  scale_range <- scale_max - scale_min
  known_sd <- is_valid_number(outcome_sd) && outcome_sd > 0
  plausible_sd <- if (is_valid_number(scale_range) && scale_range > 0) {
    scale_range / 4
  } else {
    NA_real_
  }
  sd_used <- if (known_sd) outcome_sd else plausible_sd

  if (!is_valid_number(sd_used) || sd_used <= 0) {
    return(NA_real_)
  }

  raw_effect / sd_used
}

format_effect_size <- function(x) {
  if (!is_valid_number(x)) {
    return("Enter a valid raw effect and SD or scale range.")
  }
  sprintf("%.3f SD", x)
}

custom_memory_columns <- function() {
  c(
    "simulation_id",
    "saved_at",
    "csv_file",
    "waves",
    "design_type",
    "N",
    "tau",
    "J",
    "p_comm_flooded",
    "p_indiv_given_comm",
    "tau_community",
    "tau_individual",
    "rho_y",
    "attrition_rate",
    "differential_attrition",
    "treat_prob",
    "flood_response_boost",
    "sims",
    paste0("power_", WP3_COMMUNITY_EFFECT_ORDER)
  )
}

empty_custom_memory <- function() {
  out <- as.data.frame(
    setNames(
      replicate(length(custom_memory_columns()), logical(0), simplify = FALSE),
      custom_memory_columns()
    ),
    check.names = FALSE
  )
  character_cols <- c("saved_at", "csv_file", "waves", "design_type")
  out[character_cols] <- lapply(out[character_cols], as.character)
  numeric_cols <- setdiff(
    custom_memory_columns(),
    character_cols
  )
  out[numeric_cols] <- lapply(out[numeric_cols], as.numeric)
  out
}

read_custom_memory <- function(memory_dir) {
  if (!dir.exists(memory_dir)) {
    return(empty_custom_memory())
  }

  files <- list.files(
    memory_dir,
    pattern = "^custom_power_[0-9]{8}_[0-9]{6}_sim[0-9]{3,}\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(empty_custom_memory())
  }

  rows <- lapply(files, function(path) {
    out <- tryCatch(
      read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(out)) {
      return(NULL)
    }
    missing_cols <- setdiff(custom_memory_columns(), names(out))
    for (col in missing_cols) {
      out[[col]] <- NA
    }
    out[custom_memory_columns()]
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(empty_custom_memory())
  }

  dplyr::bind_rows(rows) |>
    arrange(simulation_id, waves)
}

next_custom_simulation_id <- function(memory) {
  if (nrow(memory) == 0 || all(is.na(memory$simulation_id))) {
    return(1L)
  }
  as.integer(max(memory$simulation_id, na.rm = TRUE) + 1L)
}

custom_design_types_from_inputs <- function(waves) {
  if (is.null(waves)) {
    return(character(0))
  }

  out <- character(0)
  if ("post_only" %in% waves) {
    out <- c(out, "wp3_community_post_only")
  }
  if ("panel" %in% waves) {
    out <- c(out, "wp3_community_panel")
  }

  intersect(CUSTOM_DESIGN_TYPE_ORDER, out)
}

custom_memory_rows <- function(
  results,
  simulation_id,
  saved_at,
  csv_file,
  sims
) {
  dplyr::bind_rows(lapply(CUSTOM_DESIGN_TYPE_ORDER, function(dtype) {
    effect_order <- if (startsWith(dtype, "wp3_community")) {
      WP3_COMMUNITY_EFFECT_ORDER
    } else {
      WP3_EFFECT_ORDER
    }
    d <- results |>
      filter(
        design_type == dtype,
        outcome_scale == "Latent SD",
        term %in% effect_order
      )
    if (nrow(d) == 0) {
      return(NULL)
    }
    param <- function(col) {
      if (!col %in% names(d)) {
        return(NA_real_)
      }
      vals <- unique(stats::na.omit(d[[col]]))
      if (length(vals) == 0) {
        return(NA_real_)
      }
      vals[1]
    }

    community_powers <- setNames(
      as.list(rep(NA_real_, length(WP3_COMMUNITY_EFFECT_ORDER))),
      paste0("power_", WP3_COMMUNITY_EFFECT_ORDER)
    )
    for (effect_term in effect_order) {
      val <- d$power[d$term == effect_term]
      if (length(val) > 0) {
        target <- paste0("power_", effect_term)
        if (target %in% names(community_powers)) {
          community_powers[[target]] <- val[1]
        }
      }
    }

    as.data.frame(
      c(
        list(
          simulation_id = simulation_id,
          saved_at = saved_at,
          csv_file = csv_file,
          waves = dplyr::case_when(
            dtype == "wp3_post_only" ~ "Waves = 1",
            dtype == "wp3_panel" ~ "Waves = 2",
            dtype == "wp3_community_post_only" ~ "Community waves = 1",
            dtype == "wp3_community_panel" ~ "Community waves = 2",
            TRUE ~ dtype
          ),
          design_type = dtype,
          N = param("N"),
          tau = param("tau"),
          J = param("J"),
          p_comm_flooded = param("p_comm_flooded"),
          p_indiv_given_comm = param("p_indiv_given_comm"),
          tau_community = param("tau_community"),
          tau_individual = param("tau_individual"),
          rho_y = param("rho_y"),
          attrition_rate = param("attrition_rate"),
          differential_attrition = param("differential_attrition"),
          treat_prob = param("treat_prob"),
          flood_response_boost = param("flood_response_boost"),
          sims = sims
        ),
        community_powers
      ),
      check.names = FALSE
    )
  })) |>
    select(all_of(custom_memory_columns()))
}

# Build a compact subtitle from the sidebar filter values that are not faceted.
filter_subtitle <- function(fv, input) {
  kv <- function(cond, k, v) if (cond) paste0(k, v) else NULL
  parts <- Filter(
    Negate(is.null),
    list(
      kv(fv != "N", "N = ", input$filter_N),
      kv(fv != "rho_y", "ρ = ", input$filter_rho),
      kv(fv != "attrition_rate", "attr = ", input$filter_attrition),
      kv(
        fv != "differential_attrition",
        "survey retention penalty = ",
        input$filter_diff
      ),
      kv(fv != "flood_exposure_rate", "flood = ", input$filter_flood_rate),
      kv(
        fv != "treat_prob",
        "p(survey experiment) = ",
        input$filter_treat_prob
      ),
      kv(fv != "tau_flood", "tau_flood = ", input$filter_tau_flood),
      kv(fv != "tau_interaction", "int = ", input$filter_tau_interaction),
      kv(
        fv != "flood_response_boost",
        "flood retention boost = ",
        input$filter_flood_boost
      ),
      kv(fv != "J", "J = ", input$filter_J),
      kv(
        fv != "p_comm_flooded",
        "p_comm_flooded = ",
        input$filter_p_comm_flooded
      ),
      kv(fv != "p_indiv_given_comm", "p_indiv = ", input$filter_p_indiv),
      kv(fv != "tau_community", "tau_community = ", input$filter_tau_community),
      kv(
        fv != "tau_individual",
        "tau_individual = ",
        input$filter_tau_individual
      )
    )
  )
  if (length(parts) == 0) {
    return(NULL)
  }
  paste(parts, collapse = "  ·  ")
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Static text
# ══════════════════════════════════════════════════════════════════════════════

# ── Interpretation helpers (shown below each plot tab) ────────────────────────

interp_power <- tags$details(
  open = "open",
  tags$summary(strong(
    "How to read these plots — and what they mean for the study"
  )),
  tags$br(),
  tags$p(
    strong("What power is:"),
    " Probability of detecting a true effect
    (p < 0.05). The dashed line at 80% is a common planning benchmark."
  ),
  tags$p(
    strong("What is on the x-axis:"),
    " In the community plots, survey-effect panels use tau, community-flooding
    panels use tau_community, and individual-flooding panels use tau_individual."
  ),
  tags$p(
    strong("How to compare designs:"),
    " Post-only and panel variants are shown for the same community-level DGP.
    The primary comparison is whether the design can detect the survey effect,
    the community flooding effect, and the supplementary individual flooding
    effect."
  ),
  tags$p(
    strong("Main power drivers:"),
    " Larger N, larger absolute effect size,
    higher rho_y (for panel), and less attrition increase power. For the
    community design, J, p_comm_flooded, p_indiv_given_comm, and treat_prob also
    determine effective information."
  ),
  tags$p(
    strong("Interpretation tip:"),
    " Community flooding power depends heavily on the number of flooded
    communities; individual flooding power depends on the number of directly
    flooded respondents inside those communities."
  )
)

interp_coverage <- tags$details(
  open = "open",
  tags$summary(strong(
    "How to read these plots — and what they mean for the study"
  )),
  tags$br(),
  tags$p(
    strong("What coverage is:"),
    " Share of simulations where the 95% CI contains
    the true estimand. Target is close to 0.95."
  ),
  tags$p(
    strong("How to use this panel:"),
    " Coverage near 0.95 supports interval
    calibration under the assumed DGP. Persistent undercoverage suggests CIs are
    too narrow; persistent overcoverage suggests conservative intervals."
  ),
  tags$p(
    strong("Where risk is higher:"),
    " Coverage problems are more likely in small-N settings, high-attrition
    panel settings, and sparse community designs with few flooded communities
    or few individually flooded respondents."
  ),
  tags$p(
    strong("Practical response:"),
    " If coverage is materially off-target, consider
    adding robust/cluster-robust variants as additional design types and compare
    them directly in this app."
  )
)

interp_bias_rmse <- tags$details(
  open = "open",
  tags$summary(strong(
    "How to read these plots — and what they mean for the study"
  )),
  tags$br(),
  tags$p(
    strong("Bias:"),
    " Mean(estimate - estimand). Values near zero are expected
    under correct specification and random assignment/exposure assumptions."
  ),
  tags$p(
    strong("RMSE:"),
    " Total error metric combining variance and bias. Lower is better.
    In many scenarios, RMSE differences are driven more by variance than bias."
  ),
  tags$p(
    strong("Community-design note:"),
    " RMSE can differ sharply across the three primary effects because the
    survey effect uses the full randomized sample, the community effect is
    between communities, and the individual effect is within flooded
    communities."
  ),
  tags$p(
    strong("Attrition note:"),
    " In panel designs, differential response can induce
    estimand drift from the enrolled sample to the responder sample. This is a
    substantive design issue, not only a statistical one."
  )
)

# ── Methods section content ───────────────────────────────────────────────────

methods_content <- tagList(
  h3("Plain Language Summary"),
  tags$hr(),
  tags$p(
    "This app estimates ",
    strong("statistical power"),
    " for four design families:
    WP2 post-only, WP2 panel, WP3 post-only, and WP3 panel."
  ),
  tags$p(
    "Power is the probability of rejecting the null (p < 0.05) when a true
    effect exists. In pre-computed mode, each parameter cell is diagnosed with ",
    strong(sims_used),
    " simulations. In custom mode, the user chooses the number
    of simulations."
  ),
  tags$ul(
    tags$li(strong("WP2"), ": Survey-embedded intervention only."),
    tags$li(
      strong("WP3"),
      ": Survey intervention + natural flood exposure + interaction term."
    ),
    tags$li(strong("Post-only"), ": Uses wave-2 outcome only."),
    tags$li(
      strong("Panel"),
      ": Adds wave-1 baseline outcome (ANCOVA) and applies wave-2 attrition."
    )
  ),
  tags$p(
    "In the power plot, the x-axis depends on the selected WP3 estimand:
    survey-effect curves use ",
    tags$code("tau"),
    ", flood-effect curves use ",
    tags$code("tau_flood"),
    ", and interaction curves use ",
    tags$code("tau"),
    " with flood-effect facets. Other plot tabs continue to use ",
    tags$code("tau"),
    " as the x-axis."
  ),
  tags$br(),
  h3("Technical Description"),
  tags$hr(),
  h4("Framework"),
  tags$p(
    "The app uses ",
    tags$a("DeclareDesign", href = "https://declaredesign.org"),
    " and Monte Carlo simulation. Each design is declared as Model, Inquiry,
    Data strategy, and Answer strategy, then diagnosed repeatedly."
  ),
  h4("Model"),
  tags$p("Each simulated sample includes N respondents with:"),
  tags$ul(
    tags$li(
      "Demographics (gender, age group, education, ideology) calibrated to plausible Danish proportions."
    ),
    tags$li("Climate concern correlated with ideology."),
    tags$li(
      "Wave-1 latent outcome ",
      tags$code("Y_pre_latent"),
      " generated from demographics + noise."
    ),
    tags$li(
      "Wave-2 control potential outcome ",
      tags$code("Y0_latent = rho_y * Y_pre_latent + sqrt(1-rho_y^2) * shock"),
      " to enforce the target cross-wave correlation."
    ),
    tags$li(
      "Survey experiment effect heterogeneity ",
      tags$code("tau_z_i = tau + 0.03 * climate_high - 0.02 * ideology_right"),
      "."
    ),
    tags$li(
      "WP3 natural exposure indicator ",
      tags$code("flooded ~ Bernoulli(flood_exposure_rate)"),
      "."
    ),
    tags$li(
      "WP3 potential outcomes in a 2x2 structure over survey experiment assignment ",
      tags$code("Z"),
      " and flood exposure ",
      tags$code("F"),
      ", including flood main effect ",
      tags$code("tau_flood"),
      " and interaction effect ",
      tags$code("tau_interaction"),
      "."
    ),
    tags$li(
      "The WP3 community design assigns ",
      tags$code("N"),
      " individuals to ",
      tags$code("J"),
      " communities, draws community flooding as ",
      tags$code("Bernoulli(p_comm_flooded)"),
      ", and draws individual flooding only within flooded
      communities using ",
      tags$code("p_indiv_given_comm"),
      "."
    ),
    tags$li(
      "Community-design potential outcomes are additive in survey treatment,
      community flooding, and individual flooding. The DGP has no
      nonzero ",
      tags$code("Z x community"),
      " or ",
      tags$code("Z x individual"),
      " effect in this version; those interaction coefficients are retained as
      diagnostics rather than assigned nonzero target estimands."
    ),
    tags$li(
      "Observed Likert outcome is a coarsened 1-5 transform of the latent outcome."
    )
  ),
  h4("Inquiries (Estimands)"),
  tags$p("WP2 and WP3 track different estimands:"),
  tags$ul(
    tags$li(
      strong("WP2"),
      ": ",
      tags$code("ATE_latent"),
      " and ",
      tags$code("ATE_1_5"),
      " for survey experiment effect only."
    ),
    tags$li(
      strong("WP3"),
      ": average, subgroup, and contrast estimands for both survey experiment
      effects and flooding effects, each on latent and Likert scales."
    )
  ),
  h4("Data strategy"),
  tags$ul(
    tags$li(
      strong("Survey experiment assignment:"),
      " Complete randomization with probability ",
      tags$code("treat_prob"),
      " (not fixed at 0.50)."
    ),
    tags$li(
      strong("Outcome reveal:"),
      " The observed wave-2 outcome is selected from the 2x2 potential-outcome set based on observed ",
      tags$code("Z"),
      " and ",
      tags$code("flooded"),
      "."
    ),
    tags$li(
      strong("Panel attrition model:"),
      " In panel designs, response at wave 2 is Bernoulli with log-odds:",
      tags$code(
        "qlogis(1-attrition_rate) - differential_attrition*Z + flood_response_boost*flooded + 0.15*climate_high - 0.10*(education=='Low')"
      )
    ),
    tags$li(
      strong("Interpretation of retention parameters:"),
      " Positive ",
      tags$code("differential_attrition"),
      " lowers follow-up response among survey experiment treated respondents.
      Positive ",
      tags$code("flood_response_boost"),
      " raises follow-up response among flood-exposed respondents; negative
      values lower it."
    )
  ),
  h4("Answer strategy (Estimators)"),
  tags$ul(
    tags$li(strong("WP2 post-only:"), tags$code("Y ~ Z")),
    tags$li(strong("WP2 panel:"), tags$code("Y ~ Z + Y_pre_latent")),
    tags$li(strong("WP3 post-only:"), tags$code("Y ~ Z * flooded")),
    tags$li(strong("WP3 panel:"), tags$code("Y ~ Z * flooded + Y_pre_latent")),
    tags$li(
      strong("WP3 community post-only:"),
      tags$code(
        "lm_robust(Y ~ Z * community_flooded + Z * individual_flooded, clusters = community_id, se_type = 'CR2')"
      )
    ),
    tags$li(
      strong("WP3 community panel:"),
      tags$code(
        "lm_robust(Y ~ Z * community_flooded + Z * individual_flooded + Y_pre_latent, clusters = community_id, se_type = 'CR2')"
      )
    ),
    tags$li(
      "WP2 and standard WP3 estimators use OLS with IID standard errors. The
      community design uses community-clustered CR2 standard errors for all
      terms because IID standard errors would be invalid under the hierarchical
      exposure process."
    ),
    tags$li(
      "Effective N differs by estimand in the community design: community effects
      depend on the number of flooded communities, while individual flooding
      effects depend on individually flooded respondents inside flooded
      communities."
    ),
    tags$li(
      "The individual flooding coefficient targets the incremental flooding
      effect conditional on community flooding, ",
      tags$code("tau_individual"),
      ", not a prevalence-weighted population-average burden."
    )
  ),
  h4("Pre-computed parameter grids"),
  tags$p(
    "Default grid values are intentionally constrained to keep runtime manageable while still covering key WP3 uncertainty:"
  ),
  tags$ul(
    tags$li(
      tags$code("flood_exposure_rate"),
      " uses low-prevalence values (e.g., 0.10 and 0.20)."
    ),
    tags$li(
      tags$code("treat_prob"),
      " includes balanced and imbalanced assignment (e.g., 0.50 and 0.65)."
    ),
    tags$li(tags$code("tau_flood"), " includes at least two negative values."),
    tags$li(
      tags$code("tau_interaction"),
      " includes both negative and positive values."
    ),
    tags$li(
      tags$code("flood_response_boost"),
      " varies follow-up retention among flood-exposed respondents."
    ),
    tags$li(
      "WP3 panel attrition baseline vectors are currently fixed for runtime, but editable in ",
      tags$code("precompute.R"),
      "."
    )
  ),
  h4("Diagnosands"),
  tags$p("Computed for every (design, parameter cell) combination:"),
  tags$table(
    class = "table table-sm table-bordered",
    tags$thead(tags$tr(
      tags$th("Diagnosand"),
      tags$th("Formula"),
      tags$th("Interpretation")
    )),
    tags$tbody(
      tags$tr(
        tags$td("Power"),
        tags$td("mean(p < 0.05)"),
        tags$td("Rejection rate under the true effect")
      ),
      tags$tr(
        tags$td("Bias"),
        tags$td("mean(estimate \u2212 estimand)"),
        tags$td("Systematic over/under-estimation")
      ),
      tags$tr(
        tags$td("RMSE"),
        tags$td("\u221a mean((estimate\u2212estimand)\u00b2)"),
        tags$td("Combined accuracy (bias + variance)")
      ),
      tags$tr(
        tags$td("Coverage"),
        tags$td("mean(CI contains estimand)"),
        tags$td("95% CI nominal performance")
      ),
      tags$tr(
        tags$td("Mean SE"),
        tags$td("mean(std.error)"),
        tags$td("Average estimated uncertainty")
      ),
      tags$tr(
        tags$td("Type S"),
        tags$td("mean(p<0.05 & wrong sign)"),
        tags$td("Directional error among significant results")
      )
    )
  ),
  h4("Simulation details"),
  tags$ul(
    tags$li(
      "Simulations per cell: ",
      strong(sims_used),
      " (pre-computed grid). Custom runs use a user-specified count."
    ),
    tags$li("No bootstrap SEs on diagnosands (bootstrap_sims = FALSE)."),
    tags$li("Seed fixed at 20260413 in precompute.R for reproducibility."),
    tags$li("Pre-computed timestamp: ", precompute_ts, "."),
    tags$li(
      "Because WP3 includes interaction estimands, precision can drop sharply when flood prevalence is low or treatment split is highly imbalanced."
    )
  ),
  h4("Caveats"),
  tags$ul(
    tags$li(
      "Simulations describe operating characteristics under assumed DGPs; they do not prove identification in observed data."
    ),
    tags$li(
      "Flood exposure is currently modeled as a Bernoulli indicator with fixed prevalence, not as a geospatial process with spillovers or measurement error."
    ),
    tags$li(
      "IID OLS standard errors are used for WP2 and standard WP3. The community
      WP3 design is the clustered variant and uses CR2 standard errors at the
      community level."
    ),
    tags$li(
      "Panel estimators target effects among observed wave-2 respondents when attrition is non-ignorable."
    ),
    tags$li(
      "Power for interaction effects is usually lower than for main effects, especially when one cell in the Z x flooded table is small."
    ),
    tags$li(
      "Multiplicity adjustments are not applied; diagnosands are per-estimator."
    )
  )
)

# ── Preliminary findings content ─────────────────────────────────────────────

preliminary_findings_content <- tagList(
  h3("Preliminary Findings"),
  tags$hr(),
  tags$p(
    "These findings are based on a full parameter sweep of 9,072 design cells",
    " (1,296 WP3 community post-only + 7,776 WP3 community panel) with",
    " 500 simulations per cell (4,536,000 total simulations).",
    " All results below use the latent-SD scale.",
    " Medians are taken over all non-highlighted parameters within each grouping",
    " unless otherwise noted."
  ),
  h4("Quick simulation snapshot"),
  tags$ul(
    tags$li(
      strong("Survey effect (Z) power:"),
      " The easiest estimand to detect.",
      " Post-only median power: N = 2,000 → 0.813, N = 4,000 → 0.967,",
      " N = 6,000 → 0.996. Panel adds very little (+0.020 median gain).",
      " All tested N values provide adequate survey-effect power."
    ),
    tags$li(
      strong("Community flooding effect power:"),
      " The binding constraint for most realistic designs.",
      " Post-only crosses 0.80 only at N = 6,000 with p_comm_flooded = 0.40",
      " (median 0.882); at p_comm_flooded = 0.20, N = 6,000 post-only reaches",
      " only 0.751. Panel raises this substantially: at p_comm_flooded = 0.20,",
      " N = 4,000 panel achieves 0.695 / 0.782 / 0.921 at",
      " rho_y = 0.40 / 0.60 / 0.80.",
      " At p_comm_flooded = 0.40, N = 4,000 panel reaches 0.971 (all rho)."
    ),
    tags$li(
      strong("Individual flooding effect power:"),
      " Severely underpowered except under the most favourable parameter",
      " combination. The best cell in the grid — N = 6,000, panel,",
      " p_comm_flooded = 0.40, p_indiv_given_comm = 0.30 — reaches 0.846.",
      " At p_comm_flooded = 0.20 with p_indiv_given_comm = 0.30, N = 6,000",
      " panel achieves only 0.564. Under p_comm_flooded ≤ 0.10, power stays",
      " below 0.33 at any N. Treat this estimand as exploratory."
    ),
    tags$li(
      strong("Panel vs post-only gain (community flooding):"),
      " Overall median gain is +0.100 (Q10: +0.008, Q90: +0.356).",
      " The gain is substantial and grows sharply with rho_y",
      " (+0.044 at rho_y = 0.40, +0.130 at rho_y = 0.60, +0.284 at rho_y = 0.80)."
    ),
    tags$li(
      strong("Interaction power (Z × community flooding,",
             " Z × individual flooding):"),
      " Median power ≈ 0.05 in both post-only and panel; maximum across the",
      " full grid is 0.09. These terms are included as diagnostics; power here",
      " is at the type I error rate. Treat all interaction tests as exploratory."
    )
  ),
  h4("Key levers for one-wave vs two-wave"),
  tags$ul(
    tags$li(
      strong("rho_y is the decisive lever for the panel advantage."),
      " Median panel gain on community flooding rises from +0.044 at",
      " rho_y = 0.40 to +0.130 at rho_y = 0.60 and +0.284 at rho_y = 0.80.",
      " At rho_y ≥ 0.80, the panel design can achieve adequate community",
      " flooding power at N = 4,000 even with p_comm_flooded = 0.20."
    ),
    tags$li(
      strong("p_comm_flooded is the primary lever for community flooding power;",
             " J is nearly irrelevant."),
      " Median community flooding power shifts by only ≈0.007 across",
      " J = 300 to J = 1,200. Increasing the number of communities does not",
      " substitute for increasing the share of flooded communities."
    ),
    tags$li(
      strong("Attrition rate has negligible impact."),
      " Median community flooding power is 0.787, 0.784, and 0.782 at",
      " 10%, 30%, and 50% attrition — effectively flat."
    ),
    tags$li(
      strong("Flood retention boost (flood_response_boost) is negligible."),
      " Varying this parameter from 0 to 0.10 shifts median community",
      " flooding power by < 0.001."
    )
  ),
  h4("Design implications"),
  tags$ul(
    tags$li(
      strong("p_comm_flooded and N are the primary bottleneck levers."),
      " Post-only requires p_comm_flooded ≥ 0.40 at N = 6,000 to cross 0.80",
      " for community flooding. At p_comm_flooded = 0.20, no post-only design",
      " in the tested grid achieves 0.80."
    ),
    tags$li(
      strong("Panel design is essential when p_comm_flooded ≤ 0.20."),
      " At p_comm_flooded = 0.20 and rho_y ≥ 0.80, a panel with N = 4,000",
      " achieves 0.921. Below rho_y = 0.80 or with N = 2,000, even a panel",
      " falls short. The panel’s value is concentrated in community flooding",
      " power — it adds almost nothing for the survey effect (+0.020 gain)."
    ),
    tags$li(
      strong("Individual flooding should be designated secondary/exploratory."),
      " Only the most favourable combination (N = 6,000, panel,",
      " p_comm_flooded = 0.40, p_indiv_given_comm = 0.30) achieves 0.846.",
      " Any realistic constraint on p_comm_flooded or p_indiv_given_comm",
      " drops individual flooding power well below 0.80."
    ),
    tags$li(
      strong("Inflating J does not rescue flood-effect power."),
      " J has negligible effect on all estimands. Sampling more communities",
      " at the same p_comm_flooded yields essentially the same power.",
      " Resources are better spent increasing N or targeting higher",
      " flood-prevalence areas."
    ),
    tags$li(
      strong("Interaction tests require treating as exploratory."),
      " No design in this grid provides meaningful power for the",
      " Z × community flooding or Z × individual flooding interactions",
      " (max ≈0.09 across all cells)."
    )
  ),
  h4("Dimensions to explore next in Custom Simulation"),
  tags$ul(
    tags$li(
      "Test p_comm_flooded above 0.40 to determine whether post-only power",
      " can be achieved at N = 4,000 or lower; the tested range is not yet",
      " saturated at the top."
    ),
    tags$li(
      "Probe N between 4,000 and 6,000 at p_comm_flooded = 0.20 in the",
      " panel design (rho_y = 0.60–0.80) to locate the precise",
      " power-adequate threshold for community flooding."
    ),
    tags$li(
      "Test p_indiv_given_comm above 0.30 to map the minimum combination",
      " needed for reliable individual flooding power; this is the binding",
      " constraint for that estimand."
    ),
    tags$li(
      "If rho_y < 0.40 is plausible in the study context, test lower values",
      " (e.g., 0.20–0.30); the panel advantage may diminish substantially",
      " below the tested range."
    )
  ),
  h4("Caveat on these findings"),
  tags$p(
    "Results come from a full community-design grid sweep:",
    " N ∈ {2000, 4000, 6000},",
    " tau ∈ {0.10, 0.20},",
    " rho_y ∈ {0.40, 0.60, 0.80},",
    " J ∈ {300, 600, 1200},",
    " p_comm_flooded ∈ {0.10, 0.20, 0.40},",
    " p_indiv_given_comm ∈ {0.10, 0.30},",
    " tau_community ∈ {−0.10, −0.20},",
    " tau_individual ∈ {−0.10, −0.20},",
    " treat_prob = 0.50;",
    " panel additionally varies attrition ∈ {0.10, 0.30, 0.50} and",
    " flood_response_boost ∈ {0.00, 0.10}.",
    " 500 simulations per cell. To reproduce, run ",
    tags$code("Rscript shiny/precompute.R"),
    "."
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# 4. UI
# ══════════════════════════════════════════════════════════════════════════════

# ── Sidebar: precomputed controls ─────────────────────────────────────────────

precomputed_sidebar <- sidebar(
  width = 290,

  # ── Mode explanation ─────────────────────────────────────────────────────
  card(
    card_header("About this panel"),
    card_body(
      p(
        strong("Pre-computed mode."),
        " Results come from a grid of ",
        strong(paste(n_cells, "design cells")),
        " diagnosed before the app launched. Each cell used ",
        strong(paste(sims_used, "simulations.")),
        " Controls snap to pre-computed values — updates are instant."
      ),
      p(
        "To explore a parameter combination not in the grid, use the ",
        strong("Custom Simulation"),
        " tab. That tab runs DeclareDesign live
        and takes approximately 30\u201360 seconds."
      )
    ),
    class = "bg-light border-info"
  ),

  div(
    class = "alert alert-secondary",
    tags$strong("Scale note: "),
    "Effects are shown in latent standard-deviation units. Interpret these
    roughly like Cohen's d: 0.10 is about one tenth of a standard deviation.
    A raw Likert 1-5 difference is numerically comparable only if that Likert
    outcome has SD = 1; otherwise the raw Likert-point scale and SD scale differ."
  ),

  tags$hr(),

  # ── Display options ───────────────────────────────────────────────────────
  selectInput(
    "facet_by",
    "Facet plots by",
    choices = FACET_CHOICES,
    selected = "N"
  ),

  tags$hr(),

  # ── Parameter filters ─────────────────────────────────────────────────────
  p(em("Set values for parameters not used as the facet variable.")),

  conditionalPanel(
    condition = "input.facet_by !== 'N'",
    selectInput(
      "filter_N",
      "Sample size (N)",
      choices = N_vals,
      selected = max(N_vals)
    ),
    param_note("Number of respondents in the simulated study.")
  ),
  selectInput(
    "filter_tau",
    "Survey experiment effect (when not x-axis)",
    choices = tau_vals,
    selected = 0.10
  ),
  param_note(
    "Used to hold the survey experiment effect fixed when another effect is on the x-axis."
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'rho_y'",
    selectInput(
      "filter_rho",
      "Cross-wave correlation (\u03c1)",
      choices = rho_vals,
      selected = 0.40
    ),
    param_note(
      "Correlation between baseline and follow-up outcomes; higher values make panel adjustment more useful."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'attrition_rate'",
    selectInput(
      "filter_attrition",
      "Attrition rate (panel only)",
      choices = attr_vals,
      selected = 0.30
    ),
    param_note("Baseline probability of not completing the follow-up survey.")
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'differential_attrition'",
    selectInput(
      "filter_diff",
      "Survey experiment retention penalty (panel only)",
      choices = diff_vals,
      selected = 0.00
    ),
    param_note(
      "Log-odds penalty for follow-up response among survey experiment treated respondents; zero means no treatment-related retention difference."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'flood_exposure_rate'",
    selectInput(
      "filter_flood_rate",
      "Flood exposure rate (WP3 only)",
      choices = flood_vals,
      selected = min(flood_vals)
    ),
    param_note("Share of respondents exposed to flooding between survey waves.")
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'treat_prob'",
    selectInput(
      "filter_treat_prob",
      "Survey experiment assignment probability (WP3 only)",
      choices = treat_vals,
      selected = min(treat_vals)
    ),
    param_note(
      "Probability of assignment to the survey experiment treatment condition."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'tau_flood'",
    selectInput(
      "filter_tau_flood",
      "Flooding effect (WP3 only)",
      choices = tau_flood_vals,
      selected = min(tau_flood_vals)
    ),
    param_note(
      "Effect of flood exposure on the outcome among survey experiment controls."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'tau_interaction'",
    selectInput(
      "filter_tau_interaction",
      "Survey experiment x flooding contrast (WP3 only)",
      choices = tau_int_vals,
      selected = min(tau_int_vals)
    ),
    param_note(
      "How much the survey experiment effect differs between flood-exposed and non-exposed respondents."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'flood_response_boost'",
    selectInput(
      "filter_flood_boost",
      "Flood retention boost (WP3 panel only)",
      choices = flood_boost_vals,
      selected = if (0 %in% flood_boost_vals) 0 else min(flood_boost_vals)
    ),
    param_note(
      "Log-odds boost for follow-up response among flood-exposed respondents; negative values imply a retention penalty."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'J'",
    selectInput(
      "filter_J",
      "Number of communities (community WP3 only)",
      choices = J_vals,
      selected = if (200 %in% J_vals) 200 else max(J_vals)
    ),
    param_note(
      "Number of communities used in the hierarchical flooding design."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'p_comm_flooded'",
    selectInput(
      "filter_p_comm_flooded",
      "Flooded community share (community WP3 only)",
      choices = p_comm_flooded_vals,
      selected = min(p_comm_flooded_vals)
    ),
    param_note("Probability that a community experiences flooding.")
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'p_indiv_given_comm'",
    selectInput(
      "filter_p_indiv",
      "Individual flooding rate within flooded communities",
      choices = p_indiv_vals,
      selected = if (0.40 %in% p_indiv_vals) 0.40 else min(p_indiv_vals)
    ),
    param_note(
      "Probability of individual flooding among respondents in flooded communities."
    )
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'tau_community'",
    selectInput(
      "filter_tau_community",
      "Community flooding effect (community WP3 only)",
      choices = tau_community_vals,
      selected = if (-0.10 %in% tau_community_vals) {
        -0.10
      } else {
        min(tau_community_vals)
      }
    ),
    param_note("Community-level flood effect in latent SD units.")
  ),
  conditionalPanel(
    condition = "input.facet_by !== 'tau_individual'",
    selectInput(
      "filter_tau_individual",
      "Individual flooding effect (community WP3 only)",
      choices = tau_individual_vals,
      selected = if (-0.05 %in% tau_individual_vals) {
        -0.05
      } else {
        min(tau_individual_vals)
      }
    ),
    param_note(
      "Incremental individual flooding effect within flooded communities."
    )
  ),

  tags$hr(),
  p(tags$small(tags$em(
    "Note: WP2-only and WP3-only parameters do not apply to all designs. ",
    "Rows where a parameter is not relevant are retained unless that parameter ",
    "is used for faceting."
  )))
)

# ── Main tabs: pre-computed dashboard ─────────────────────────────────────────

dashboard_main <- navset_tab(
  id = "dashboard_tabs",

  nav_panel(
    "Power",
    div(
      class = "alert alert-secondary",
      tags$strong("Effective sample-size diagnostic: "),
      "The community effect is identified between communities, so effective N
      depends on the number of flooded communities. The individual flooding
      effect is identified within flooded communities. This is why the two
      flooding outputs can have very different precision under the same N."
    ),
    wp3_community_effect_grid_ui("community_power"),
    interp_power
  ),

  nav_panel(
    "Diagnostics",
    div(
      class = "alert alert-secondary",
      tags$strong("Diagnostics are secondary. "),
      "Use these mainly as checks on interval calibration and estimator behavior;
      power is the primary planning metric."
    ),
    accordion(
      accordion_panel(
        "Coverage",
        wp3_community_effect_grid_ui("community_coverage", height = "240px"),
        interp_coverage
      ),
      accordion_panel(
        "Bias",
        wp3_community_effect_grid_ui("community_bias", height = "240px")
      ),
      accordion_panel(
        "RMSE",
        wp3_community_effect_grid_ui("community_rmse", height = "240px"),
        interp_bias_rmse
      ),
      open = FALSE
    )
  )
)

# ── Custom simulation panel ───────────────────────────────────────────────────

custom_panel <- fluidRow(
  column(
    4,
    card(
      card_header("Custom parameters"),
      card_body(
        div(
          class = "alert alert-warning",
          tags$strong("Heads up:"),
          " This runs DeclareDesign live in your R
          session. Expect ",
          tags$strong("30\u201360 seconds"),
          " per design type.
          The app will show a progress bar and remain unresponsive until
          complete. This is normal."
        ),
        p(
          "Set any combination of parameters — including values not in the
          pre-computed grid. Results will appear on the right as plots and
          a numerical table."
        ),
        div(
          class = "alert alert-info",
          tags$strong("Suggested values to test:"),
          tags$ul(
            tags$li(
              "Sample size: ",
              tags$code("4000"),
              ", ",
              tags$code("6000"),
              ", ",
              tags$code("8000")
            ),
            tags$li(
              "Survey experiment effect: ",
              tags$code("0.05"),
              ", ",
              tags$code("0.10"),
              ", ",
              tags$code("0.15"),
              " SD"
            ),
            tags$li(
              "Flooding effect: ",
              tags$code("-0.10"),
              ", ",
              tags$code("-0.05"),
              ", ",
              tags$code("0.00"),
              " SD"
            ),
            tags$li(
              "Survey experiment x flooding contrast: ",
              tags$code("-0.08"),
              ", ",
              tags$code("-0.04"),
              ", ",
              tags$code("0.04"),
              ", ",
              tags$code("0.08"),
              " SD"
            ),
            tags$li(
              "Flood exposure rate: ",
              tags$code("0.05"),
              ", ",
              tags$code("0.10"),
              ", ",
              tags$code("0.20")
            ),
            tags$li(
              "Retention assumptions: attrition ",
              tags$code("0.20"),
              "-",
              tags$code("0.40"),
              ", survey experiment retention penalty ",
              tags$code("0.00"),
              ", flood retention boost ",
              tags$code("-0.10"),
              ", ",
              tags$code("0.00"),
              ", ",
              tags$code("0.10")
            ),
            tags$li(
              "Use ",
              tags$code("200"),
              " simulations for quick screening and ",
              tags$code("1000"),
              "+ before treating a scenario as stable."
            )
          )
        ),
        tags$hr(),
        numericInput(
          "custom_N",
          "Sample size (N)",
          value = 4000,
          min = 100,
          max = 20000,
          step = 100
        ),
        param_note("Number of respondents in the simulated study."),
        accordion(
          accordion_panel(
            "Effect-size helper",
            p(
              "Convert a raw outcome difference into the standardized effects
              used by the simulation. If the outcome SD is unknown, the helper
              uses range / 4 as a common planning approximation."
            ),
            numericInput(
              "helper_raw_effect",
              "Raw effect",
              value = 0.10,
              min = -100,
              max = 100,
              step = 0.01
            ),
            param_note(
              "Difference on the original outcome scale, such as Likert points."
            ),
            numericInput(
              "helper_outcome_sd",
              "Known outcome SD",
              value = NA,
              min = 0,
              max = 100,
              step = 0.01
            ),
            param_note("Leave blank if the SD is unknown."),
            fluidRow(
              column(
                6,
                numericInput(
                  "helper_scale_min",
                  "Scale minimum",
                  value = 1,
                  min = -100,
                  max = 100,
                  step = 1
                )
              ),
              column(
                6,
                numericInput(
                  "helper_scale_max",
                  "Scale maximum",
                  value = 5,
                  min = -100,
                  max = 100,
                  step = 1
                )
              )
            ),
            uiOutput("effect_size_helper_summary"),
            fluidRow(
              column(
                4,
                actionButton(
                  "use_helper_tau",
                  "Survey",
                  icon = icon("arrow-right"),
                  class = "btn-outline-primary btn-sm w-100"
                )
              ),
              column(
                4,
                actionButton(
                  "use_helper_tau_community",
                  "Community",
                  icon = icon("arrow-right"),
                  class = "btn-outline-primary btn-sm w-100"
                )
              ),
              column(
                4,
                actionButton(
                  "use_helper_tau_individual",
                  "Individual",
                  icon = icon("arrow-right"),
                  class = "btn-outline-primary btn-sm w-100"
                )
              )
            ),
            param_note(
              "Buttons copy the standardized value into the matching effect input below."
            ),
            value = "effect-size-helper"
          ),
          open = FALSE
        ),
        numericInput(
          "custom_tau",
          "Survey experiment effect (SD units)",
          value = 0.10,
          min = 0,
          max = 2,
          step = 0.01
        ),
        param_note(
          "Average effect of the embedded survey experiment, in latent SD units."
        ),
        numericInput(
          "custom_rho",
          "Cross-wave correlation (\u03c1)",
          value = 0.40,
          min = -1,
          max = 1,
          step = 0.05
        ),
        param_note(
          "Correlation between baseline and follow-up outcomes; higher values make panel adjustment more useful."
        ),
        numericInput(
          "custom_attr",
          "Attrition rate",
          value = 0.30,
          min = 0,
          max = 0.9,
          step = 0.05
        ),
        param_note(
          "Baseline probability of not completing the follow-up survey."
        ),
        numericInput(
          "custom_diff",
          "Survey experiment retention penalty (log-odds)",
          value = 0.00,
          min = 0,
          max = 1,
          step = 0.05
        ),
        param_note(
          "Positive values lower follow-up response among survey experiment treated respondents; zero means no treatment-related retention difference."
        ),
        tags$hr(),
        h4("Community flooding inputs"),
        numericInput(
          "custom_J",
          "Number of communities",
          value = 200,
          min = 10,
          max = 5000,
          step = 10
        ),
        param_note(
          "Number of communities in the hierarchical flooding design."
        ),
        numericInput(
          "custom_p_comm_flooded",
          "Share of flooded communities",
          value = 0.15,
          min = 0,
          max = 1,
          step = 0.01
        ),
        param_note("Probability that a community experiences flooding."),
        numericInput(
          "custom_p_indiv",
          "Individual flooding rate within flooded communities",
          value = 0.40,
          min = 0,
          max = 1,
          step = 0.01
        ),
        param_note(
          "Probability of individual flooding among respondents in flooded communities."
        ),
        numericInput(
          "custom_tau_community",
          "Community flooding effect (SD units)",
          value = -0.10,
          min = -2,
          max = 2,
          step = 0.01
        ),
        param_note("Community-level flood effect in latent SD units."),
        numericInput(
          "custom_tau_individual",
          "Individual flooding effect (SD units)",
          value = -0.05,
          min = -2,
          max = 2,
          step = 0.01
        ),
        param_note(
          "Supplementary individual flooding effect within flooded communities."
        ),
        numericInput(
          "custom_treat_prob",
          "Survey experiment assignment probability",
          value = 0.50,
          min = 0.05,
          max = 0.95,
          step = 0.01
        ),
        param_note(
          "Probability of assignment to the survey experiment treatment condition."
        ),
        numericInput(
          "custom_flood_boost",
          "Flood retention boost (log-odds)",
          value = 0.00,
          min = -1,
          max = 1,
          step = 0.01
        ),
        param_note(
          "Positive values raise follow-up response among flood-exposed respondents; negative values lower it."
        ),
        numericInput(
          "custom_sims",
          "Simulations",
          value = 200,
          min = 50,
          max = 2000,
          step = 50
        ),
        param_note("Monte Carlo repetitions for each selected design."),
        checkboxGroupInput(
          "custom_designs",
          "Waves",
          choices = CUSTOM_WAVE_CHOICES,
          selected = unname(CUSTOM_WAVE_CHOICES)
        ),
        param_note(
          "Choose whether to diagnose a one-wave post-only design, a two-wave pre-post design, or both."
        ),
        actionButton(
          "run_custom",
          "Run simulation",
          class = "btn-primary btn-lg w-100",
          icon = icon("play")
        )
      )
    )
  ),
  column(
    8,
    uiOutput("custom_status"),
    conditionalPanel(
      condition = "output.custom_has_results",
      navset_tab(
        nav_panel(
          "Power",
          h4(
            "Community design: survey, community flooding, individual flooding"
          ),
          wp3_community_effect_grid_ui("custom_community_power")
        ),
        nav_panel(
          "Diagnostics",
          div(
            class = "alert alert-secondary",
            tags$strong("Diagnostics are secondary. "),
            "Use these mainly as checks on interval calibration and estimator behavior."
          ),
          accordion(
            accordion_panel(
              "Coverage",
              h4(
                "Community design: survey, community flooding, individual flooding"
              ),
              wp3_community_effect_grid_ui(
                "custom_community_coverage",
                height = "240px"
              )
            ),
            accordion_panel(
              "Bias",
              h4(
                "Community design: survey, community flooding, individual flooding"
              ),
              wp3_community_effect_grid_ui(
                "custom_community_bias",
                height = "240px"
              )
            ),
            accordion_panel(
              "RMSE",
              h4(
                "Community design: survey, community flooding, individual flooding"
              ),
              wp3_community_effect_grid_ui(
                "custom_community_rmse",
                height = "240px"
              )
            ),
            open = FALSE
          )
        ),
        nav_panel("Numerical results", DTOutput("custom_table"))
      )
    ),
    tags$hr(),
    h4("Simulation memory"),
    uiOutput("custom_memory_status"),
    DTOutput("custom_memory_table")
  )
)

# ── Full UI ───────────────────────────────────────────────────────────────────

ui <- page_navbar(
  title = "WP2/WP3 Power Analysis",
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_google("Source Sans Pro")
  ),
  fillable = FALSE,

  nav_panel(
    "Power Dashboard",
    layout_sidebar(
      sidebar = precomputed_sidebar,
      dashboard_main
    )
  ),

  nav_panel(
    "Custom Simulation",
    div(class = "container-fluid mt-3", custom_panel)
  ),

  nav_panel(
    "Methods",
    div(class = "container mt-3", methods_content)
  ),

  nav_panel(
    "Preliminary findings",
    div(class = "container mt-3", preliminary_findings_content)
  ),

  nav_panel(
    "Interpretation Guide",
    div(
      class = "container mt-3",
      h3("Interpretation Guide"),
      p(
        "Each section below explains one diagnosand: what it measures, how to
        read the plot, and what it implies for study design decisions."
      ),
      tags$hr(),
      h4("Power"),
      interp_power,
      tags$br(),
      h4("Coverage"),
      interp_coverage,
      tags$br(),
      h4("Bias & RMSE"),
      interp_bias_rmse
    )
  ),

  nav_spacer(),
  nav_item(tags$small(tags$em(
    paste0(
      "Grid: ",
      length(N_vals),
      " N \u00d7 ",
      length(tau_vals),
      " \u03c4 \u00d7 ",
      length(rho_vals),
      " \u03c1 \u00d7 ",
      length(flood_vals),
      " flood-rate \u00d7 ",
      length(treat_vals),
      " treat-prob \u00d7 ",
      length(tau_flood_vals),
      " flood-effect \u00d7 ",
      length(tau_int_vals),
      " interaction \u00d7 ",
      length(flood_boost_vals),
      " flood-retention boost | ",
      sims_used,
      " sims | ",
      precompute_ts
    )
  )))
)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Server
# ══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {
  custom_memory_dir <- "custom_simulation_memory"
  if (!dir.exists(custom_memory_dir)) {
    dir.create(custom_memory_dir, recursive = TRUE, showWarnings = FALSE)
  }
  custom_memory <- reactiveVal(read_custom_memory(custom_memory_dir))

  helper_effect_sd <- reactive({
    effect_size_helper_value(
      raw_effect = input$helper_raw_effect,
      outcome_sd = input$helper_outcome_sd,
      scale_min = input$helper_scale_min,
      scale_max = input$helper_scale_max
    )
  })

  output$effect_size_helper_summary <- renderUI({
    raw_effect <- input$helper_raw_effect
    outcome_sd <- input$helper_outcome_sd
    scale_min <- input$helper_scale_min
    scale_max <- input$helper_scale_max
    scale_range <- scale_max - scale_min
    known_sd <- is_valid_number(outcome_sd) && outcome_sd > 0

    if (!is_valid_number(raw_effect)) {
      return(div(class = "alert alert-warning", "Enter a valid raw effect."))
    }

    if (known_sd) {
      return(div(
        class = "alert alert-info",
        tags$strong("Standardized effect: "),
        format_effect_size(helper_effect_sd()),
        tags$br(),
        tags$small("Calculated as raw effect / known outcome SD.")
      ))
    }

    if (!is_valid_number(scale_range) || scale_range <= 0) {
      return(div(
        class = "alert alert-warning",
        "Enter a valid scale range or known outcome SD."
      ))
    }

    common_sd <- scale_range / 4
    uniform_sd <- scale_range / sqrt(12)
    max_sd <- scale_range / 2

    div(
      class = "alert alert-info",
      tags$strong("Planning estimate: "),
      format_effect_size(raw_effect / common_sd),
      tags$br(),
      tags$small(
        "Using common SD approx. range / 4 = ",
        sprintf("%.3f", common_sd),
        ". Uniform-spread SD would be ",
        sprintf("%.3f", uniform_sd),
        " (",
        format_effect_size(raw_effect / uniform_sd),
        "); theoretical maximum SD is ",
        sprintf("%.3f", max_sd),
        " (",
        format_effect_size(raw_effect / max_sd),
        ")."
      )
    )
  })

  observeEvent(input$use_helper_tau, {
    val <- helper_effect_sd()
    req(is.finite(val))
    updateNumericInput(session, "custom_tau", value = round(abs(val), 3))
  })

  observeEvent(input$use_helper_tau_community, {
    val <- helper_effect_sd()
    req(is.finite(val))
    updateNumericInput(session, "custom_tau_community", value = round(val, 3))
  })

  observeEvent(input$use_helper_tau_individual, {
    val <- helper_effect_sd()
    req(is.finite(val))
    updateNumericInput(session, "custom_tau_individual", value = round(val, 3))
  })

  # ── Reactive: filter pre-computed data ─────────────────────────────────────

  plot_data <- reactive({
    d <- results

    # Scale filter
    d <- filter(d, outcome_scale == "Latent SD")

    fv <- input$facet_by

    # If faceting by a parameter not used by some designs, those rows are NA
    # and cannot be faceted meaningfully.
    if (
      fv %in%
        c(
          "attrition_rate",
          "differential_attrition",
          "flood_exposure_rate",
          "treat_prob",
          "tau_flood",
          "tau_interaction",
          "flood_response_boost",
          "J",
          "p_comm_flooded",
          "p_indiv_given_comm",
          "tau_community",
          "tau_individual"
        )
    ) {
      d <- filter(d, !is.na(.data[[fv]]))
    }

    # Filter non-facet parameters to the selected value
    if (fv != "N") {
      d <- filter(d, N == as.numeric(input$filter_N))
    }
    if (fv != "rho_y") {
      d <- filter(d, rho_y == as.numeric(input$filter_rho))
    }
    if (fv != "attrition_rate") {
      d <- filter(
        d,
        is.na(attrition_rate) |
          attrition_rate == as.numeric(input$filter_attrition)
      )
    }
    if (fv != "differential_attrition") {
      d <- filter(
        d,
        is.na(differential_attrition) |
          differential_attrition == as.numeric(input$filter_diff)
      )
    }
    if (fv != "flood_exposure_rate") {
      d <- filter(
        d,
        is.na(flood_exposure_rate) |
          flood_exposure_rate == as.numeric(input$filter_flood_rate)
      )
    }
    if (fv != "treat_prob") {
      d <- filter(
        d,
        is.na(treat_prob) |
          treat_prob == as.numeric(input$filter_treat_prob)
      )
    }
    if (fv != "tau_flood") {
      d <- filter(
        d,
        is.na(tau_flood) |
          tau_flood == as.numeric(input$filter_tau_flood)
      )
    }
    if (fv != "tau_interaction") {
      d <- filter(
        d,
        is.na(tau_interaction) |
          tau_interaction == as.numeric(input$filter_tau_interaction)
      )
    }
    if (fv != "flood_response_boost") {
      d <- filter(
        d,
        is.na(flood_response_boost) |
          flood_response_boost == as.numeric(input$filter_flood_boost)
      )
    }
    if (fv != "J") {
      d <- filter(
        d,
        is.na(J) |
          J == as.numeric(input$filter_J)
      )
    }
    if (fv != "p_comm_flooded") {
      d <- filter(
        d,
        is.na(p_comm_flooded) |
          p_comm_flooded == as.numeric(input$filter_p_comm_flooded)
      )
    }
    if (fv != "p_indiv_given_comm") {
      d <- filter(
        d,
        is.na(p_indiv_given_comm) |
          p_indiv_given_comm == as.numeric(input$filter_p_indiv)
      )
    }
    if (fv != "tau_community") {
      d <- filter(
        d,
        is.na(tau_community) |
          tau_community == as.numeric(input$filter_tau_community)
      )
    }
    if (fv != "tau_individual") {
      d <- filter(
        d,
        is.na(tau_individual) |
          tau_individual == as.numeric(input$filter_tau_individual)
      )
    }

    d
  })

  # ── Shared plot helper ──────────────────────────────────────────────────────

  add_facet <- function(p, fv) {
    if (fv == "none") {
      return(p)
    }
    label_fn <- function(value) paste0(FACET_LABELS[fv], " = ", value)
    p + facet_wrap(as.formula(paste("~", fv)), labeller = as_labeller(label_fn))
  }

  base_theme <- function() {
    theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }

  # ── WP3 effect plots ────────────────────────────────────────────────────────

  wp3_facet_vars <- c(
    "N",
    "rho_y",
    "attrition_rate",
    "differential_attrition",
    "flood_exposure_rate",
    "treat_prob",
    "tau_flood",
    "tau_interaction",
    "flood_response_boost"
  )

  wp3_effect_data <- function(effect_term) {
    x_var <- unname(WP3_EFFECT_X[[effect_term]])
    d <- results |>
      filter(
        design_type %in% c("wp3_post_only", "wp3_panel"),
        outcome_scale == "Latent SD",
        term == effect_term
      ) |>
      mutate(
        design_variant = ifelse(
          design_type == "wp3_post_only",
          "Post-only",
          "Panel"
        )
      )

    fv <- input$facet_by

    if (fv != "N") {
      d <- filter(d, N == as.numeric(input$filter_N))
    }
    if (fv != "rho_y") {
      d <- filter(d, rho_y == as.numeric(input$filter_rho))
    }
    if (fv != "attrition_rate") {
      d <- filter(
        d,
        is.na(attrition_rate) |
          attrition_rate == as.numeric(input$filter_attrition)
      )
    }
    if (fv != "differential_attrition") {
      d <- filter(
        d,
        is.na(differential_attrition) |
          differential_attrition == as.numeric(input$filter_diff)
      )
    }
    if (fv != "flood_exposure_rate") {
      d <- filter(d, flood_exposure_rate == as.numeric(input$filter_flood_rate))
    }
    if (fv != "treat_prob") {
      d <- filter(d, treat_prob == as.numeric(input$filter_treat_prob))
    }
    if (x_var != "tau") {
      d <- filter(d, tau == as.numeric(input$filter_tau))
    }
    if (x_var != "tau_flood" && fv != "tau_flood") {
      d <- filter(d, tau_flood == as.numeric(input$filter_tau_flood))
    }
    if (x_var != "tau_interaction" && fv != "tau_interaction") {
      d <- filter(
        d,
        tau_interaction == as.numeric(input$filter_tau_interaction)
      )
    }
    if (fv != "flood_response_boost") {
      d <- filter(
        d,
        is.na(flood_response_boost) |
          flood_response_boost == as.numeric(input$filter_flood_boost)
      )
    }

    d
  }

  wp3_effect_plot <- function(
    effect_term,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    x_var <- unname(WP3_EFFECT_X[[effect_term]])
    fv <- input$facet_by
    d <- wp3_effect_data(effect_term)
    req(nrow(d) > 0)

    p <- ggplot(
      d,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = design_variant,
        linetype = design_variant,
        shape = design_variant,
        group = design_variant
      )
    )

    if (!is.null(hline)) {
      p <- p +
        geom_hline(
          yintercept = hline,
          linetype = "dashed",
          color = "grey40",
          linewidth = 0.7
        ) +
        annotate(
          "text",
          x = min(d[[x_var]], na.rm = TRUE),
          y = hline + 0.02,
          label = hline_label,
          hjust = 0,
          size = 3.2,
          color = "grey40"
        )
    }

    p <- p +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.5) +
      scale_color_manual(
        values = c("Post-only" = "#1f77b4", "Panel" = "#d62728"),
        name = "WP3 design"
      ) +
      scale_linetype_manual(
        values = c("Post-only" = "solid", "Panel" = "dashed"),
        name = "WP3 design"
      ) +
      scale_shape_manual(
        values = c("Post-only" = 16L, "Panel" = 17L),
        name = "WP3 design"
      ) +
      labs(
        x = WP3_EFFECT_X_LABELS[[x_var]],
        y = y_label,
        title = paste(title_prefix, WP3_EFFECT_LABELS[[effect_term]]),
        subtitle = filter_subtitle(fv, input)
      ) +
      base_theme()

    if (!is.null(y_limits) || percent_y) {
      p <- p +
        scale_y_continuous(
          labels = if (percent_y) percent_format(accuracy = 1) else waiver()
        )
    }
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }

    if (fv != x_var && fv %in% wp3_facet_vars) {
      add_facet(p, fv)
    } else {
      p
    }
  }

  register_wp3_effect_plots <- function(
    prefix,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    for (effect_term in WP3_EFFECT_ORDER) {
      local({
        effect_local <- effect_term
        output[[paste0(prefix, "_", effect_local)]] <- renderPlot({
          wp3_effect_plot(
            effect_term = effect_local,
            y_var = y_var,
            y_label = y_label,
            title_prefix = title_prefix,
            hline = hline,
            hline_label = hline_label,
            y_limits = y_limits,
            percent_y = percent_y
          )
        })
      })
    }
  }

  register_wp3_effect_plots(
    "power",
    "power",
    "Power",
    "Power:",
    hline = 0.80,
    hline_label = "80% target",
    y_limits = c(0, 1),
    percent_y = TRUE
  )

  register_wp3_effect_plots(
    "coverage",
    "coverage",
    "95% CI coverage",
    "Coverage:",
    hline = 0.95,
    hline_label = "95%",
    y_limits = c(0.80, 1),
    percent_y = TRUE
  )

  register_wp3_effect_plots(
    "bias",
    "bias",
    "Bias (estimate - estimand)",
    "Bias:",
    hline = 0
  )

  register_wp3_effect_plots(
    "rmse",
    "rmse",
    "RMSE",
    "RMSE:"
  )

  # ── WP3 community effect plots ─────────────────────────────────────────────

  community_facet_vars <- c(
    "N",
    "rho_y",
    "attrition_rate",
    "differential_attrition",
    "treat_prob",
    "flood_response_boost",
    "J",
    "p_comm_flooded",
    "p_indiv_given_comm",
    "tau_community",
    "tau_individual"
  )

  wp3_community_effect_data <- function(effect_term) {
    x_var <- unname(WP3_COMMUNITY_EFFECT_X[[effect_term]])
    d <- results |>
      filter(
        design_type %in%
          c("wp3_community_post_only", "wp3_community_panel"),
        outcome_scale == "Latent SD",
        term == effect_term
      ) |>
      mutate(
        design_variant = ifelse(
          design_type == "wp3_community_post_only",
          "Post-only",
          "Panel"
        )
      )

    fv <- input$facet_by

    if (fv != "N") {
      d <- filter(d, N == as.numeric(input$filter_N))
    }
    if (fv != "rho_y") {
      d <- filter(d, rho_y == as.numeric(input$filter_rho))
    }
    if (fv != "attrition_rate") {
      d <- filter(
        d,
        is.na(attrition_rate) |
          attrition_rate == as.numeric(input$filter_attrition)
      )
    }
    if (fv != "differential_attrition") {
      d <- filter(
        d,
        is.na(differential_attrition) |
          differential_attrition == as.numeric(input$filter_diff)
      )
    }
    if (fv != "treat_prob") {
      d <- filter(d, treat_prob == as.numeric(input$filter_treat_prob))
    }
    if (fv != "flood_response_boost") {
      d <- filter(
        d,
        is.na(flood_response_boost) |
          flood_response_boost == as.numeric(input$filter_flood_boost)
      )
    }
    if (fv != "J") {
      d <- filter(d, J == as.numeric(input$filter_J))
    }
    if (fv != "p_comm_flooded") {
      d <- filter(
        d,
        p_comm_flooded == as.numeric(input$filter_p_comm_flooded)
      )
    }
    if (fv != "p_indiv_given_comm") {
      d <- filter(
        d,
        p_indiv_given_comm == as.numeric(input$filter_p_indiv)
      )
    }
    if (x_var != "tau") {
      d <- filter(d, tau == as.numeric(input$filter_tau))
    }
    if (x_var != "tau_community" && fv != "tau_community") {
      d <- filter(
        d,
        tau_community == as.numeric(input$filter_tau_community)
      )
    }
    if (x_var != "tau_individual" && fv != "tau_individual") {
      d <- filter(
        d,
        tau_individual == as.numeric(input$filter_tau_individual)
      )
    }

    d
  }

  wp3_community_effect_plot <- function(
    effect_term,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    x_var <- unname(WP3_COMMUNITY_EFFECT_X[[effect_term]])
    fv <- input$facet_by
    d <- wp3_community_effect_data(effect_term)
    req(nrow(d) > 0)

    p <- ggplot(
      d,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = design_variant,
        linetype = design_variant,
        shape = design_variant,
        group = design_variant
      )
    )

    if (!is.null(hline)) {
      p <- p +
        geom_hline(
          yintercept = hline,
          linetype = "dashed",
          color = "grey40",
          linewidth = 0.7
        ) +
        annotate(
          "text",
          x = min(d[[x_var]], na.rm = TRUE),
          y = hline + 0.02,
          label = hline_label,
          hjust = 0,
          size = 3.2,
          color = "grey40"
        )
    }

    p <- p +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.5) +
      scale_color_manual(
        values = c("Post-only" = "#1f77b4", "Panel" = "#d62728"),
        name = "WP3 community design"
      ) +
      scale_linetype_manual(
        values = c("Post-only" = "solid", "Panel" = "dashed"),
        name = "WP3 community design"
      ) +
      scale_shape_manual(
        values = c("Post-only" = 16L, "Panel" = 17L),
        name = "WP3 community design"
      ) +
      labs(
        x = WP3_COMMUNITY_EFFECT_X_LABELS[[x_var]],
        y = y_label,
        title = paste(
          title_prefix,
          WP3_COMMUNITY_EFFECT_LABELS[[effect_term]]
        ),
        subtitle = filter_subtitle(fv, input)
      ) +
      base_theme()

    if (!is.null(y_limits) || percent_y) {
      p <- p +
        scale_y_continuous(
          labels = if (percent_y) percent_format(accuracy = 1) else waiver()
        )
    }
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }

    if (fv != x_var && fv %in% community_facet_vars) {
      add_facet(p, fv)
    } else {
      p
    }
  }

  register_wp3_community_effect_plots <- function(
    prefix,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    for (effect_term in WP3_COMMUNITY_EFFECT_ORDER) {
      local({
        effect_local <- effect_term
        output[[effect_output_id(prefix, effect_local)]] <- renderPlot({
          wp3_community_effect_plot(
            effect_term = effect_local,
            y_var = y_var,
            y_label = y_label,
            title_prefix = title_prefix,
            hline = hline,
            hline_label = hline_label,
            y_limits = y_limits,
            percent_y = percent_y
          )
        })
      })
    }
  }

  register_wp3_community_effect_plots(
    "community_power",
    "power",
    "Power",
    "Power:",
    hline = 0.80,
    hline_label = "80% target",
    y_limits = c(0, 1),
    percent_y = TRUE
  )
  register_wp3_community_effect_plots(
    "community_coverage",
    "coverage",
    "95% CI coverage",
    "Coverage:",
    hline = 0.95,
    hline_label = "95%",
    y_limits = c(0.80, 1),
    percent_y = TRUE
  )
  register_wp3_community_effect_plots(
    "community_bias",
    "bias",
    "Bias (estimate - estimand)",
    "Bias:",
    hline = 0
  )
  register_wp3_community_effect_plots(
    "community_rmse",
    "rmse",
    "RMSE",
    "RMSE:"
  )

  # ── Custom simulation ─────────────────────────────────────────────────────

  custom_results <- eventReactive(input$run_custom, {
    req(length(input$custom_designs) > 0)

    selected_design_types <- custom_design_types_from_inputs(
      input$custom_designs
    )
    req(length(selected_design_types) > 0)

    n_designs <- length(selected_design_types)

    withProgress(
      message = paste0(
        "Running DeclareDesign (",
        n_designs,
        " design",
        if (n_designs > 1) "s" else "",
        ")..."
      ),
      value = 0,
      {
        results_list <- lapply(seq_along(selected_design_types), function(i) {
          dtype <- selected_design_types[[i]]
          dlabel <- DESIGN_LABEL_LOOKUP[[dtype]]
          incProgress(
            amount = 0.1,
            detail = paste("Building:", dlabel)
          )
          d <- build_wp3_community_design(
            N = input$custom_N,
            J = input$custom_J,
            p_comm_flooded = input$custom_p_comm_flooded,
            p_indiv_given_comm = input$custom_p_indiv,
            tau = input$custom_tau,
            tau_community = input$custom_tau_community,
            tau_individual = input$custom_tau_individual,
            rho_y = input$custom_rho,
            attrition_rate = input$custom_attr,
            differential_attrition = input$custom_diff,
            treat_prob = input$custom_treat_prob,
            flood_response_boost = input$custom_flood_boost,
            type = dtype
          )
          incProgress(
            amount = 0.35,
            detail = paste(
              "Diagnosing:",
              dlabel,
              "(",
              input$custom_sims,
              "sims)"
            )
          )
          diag <- diagnose_design(
            d,
            diagnosands = diagnosands_power,
            sims = input$custom_sims,
            bootstrap_sims = FALSE
          )
          incProgress(0.05, detail = "Collecting results...")
          diag$diagnosands_df |>
            mutate(
              N = input$custom_N,
              tau = input$custom_tau,
              rho_y = input$custom_rho,
              attrition_rate = input$custom_attr,
              differential_attrition = input$custom_diff,
              flood_exposure_rate = NA_real_,
              treat_prob = input$custom_treat_prob,
              tau_flood = NA_real_,
              tau_interaction = NA_real_,
              J = input$custom_J,
              p_comm_flooded = input$custom_p_comm_flooded,
              p_indiv_given_comm = input$custom_p_indiv,
              tau_community = input$custom_tau_community,
              tau_individual = input$custom_tau_individual,
              flood_response_boost = input$custom_flood_boost,
              design_type = dtype,
              design_label_display = dlabel,
              estimator_label = dplyr::recode(estimator, !!!ESTIMATOR_LABELS),
              outcome_scale = ifelse(
                grepl("1_5", estimator),
                "Likert 1-5",
                "Latent SD"
              )
            )
        })
        run_results <- bind_rows(results_list)
        simulation_id <- next_custom_simulation_id(custom_memory())
        saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        file_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        csv_file <- sprintf(
          "custom_power_%s_sim%03d.csv",
          file_stamp,
          simulation_id
        )
        memory_rows <- custom_memory_rows(
          results = run_results,
          simulation_id = simulation_id,
          saved_at = saved_at,
          csv_file = csv_file,
          sims = input$custom_sims
        )
        if (nrow(memory_rows) > 0) {
          write.csv(
            memory_rows,
            file.path(custom_memory_dir, csv_file),
            row.names = FALSE
          )
          custom_memory(bind_rows(custom_memory(), memory_rows))
        }

        attr(run_results, "simulation_id") <- simulation_id
        attr(run_results, "csv_file") <- csv_file
        run_results
      }
    )
  })

  output$custom_has_results <- reactive({
    !is.null(custom_results()) && nrow(custom_results()) > 0
  })
  outputOptions(output, "custom_has_results", suspendWhenHidden = FALSE)

  output$custom_status <- renderUI({
    if (is.null(custom_results())) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        " Set parameters and click ",
        strong("Run simulation"),
        " to generate results."
      )
    } else {
      r <- custom_results()
      simulation_id <- attr(r, "simulation_id")
      csv_file <- attr(r, "csv_file")
      first_non_na <- function(col) {
        vals <- unique(stats::na.omit(r[[col]]))
        if (length(vals) == 0) {
          return(NA_real_)
        }
        vals[1]
      }
      params <- sprintf(
        paste0(
          "N = %d, tau_survey = %.2f, ",
          "rho = %.2f, attrition = %.0f%%, diff. attrition = %.2f, ",
          "survey_assignment_prob = %.2f, ",
          "J = %.0f, p_comm_flooded = %.2f, p_indiv = %.2f, ",
          "tau_community = %.2f, tau_individual = %.2f, ",
          "flood_retention_boost = %.2f, sims = %d"
        ),
        first_non_na("N"),
        first_non_na("tau"),
        first_non_na("rho_y"),
        100 * first_non_na("attrition_rate"),
        first_non_na("differential_attrition"),
        first_non_na("treat_prob"),
        first_non_na("J"),
        first_non_na("p_comm_flooded"),
        first_non_na("p_indiv_given_comm"),
        first_non_na("tau_community"),
        first_non_na("tau_individual"),
        first_non_na("flood_response_boost"),
        input$custom_sims
      )
      div(
        class = "alert alert-success",
        icon("check-circle"),
        " Simulation ",
        simulation_id,
        " complete",
        if (!is.null(csv_file)) paste0(" and saved to ", csv_file) else "",
        ". Parameters: ",
        em(params)
      )
    }
  })

  custom_wp3_community_effect_data <- function(effect_term) {
    r <- custom_results()
    req(!is.null(r))

    r |>
      filter(
        outcome_scale == "Latent SD",
        design_type %in%
          c("wp3_community_post_only", "wp3_community_panel"),
        term == effect_term
      ) |>
      mutate(
        design_variant = dplyr::recode(
          design_type,
          wp3_community_post_only = "Post-only",
          wp3_community_panel = "Panel"
        )
      )
  }

  custom_wp3_community_effect_plot <- function(
    effect_term,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    d <- custom_wp3_community_effect_data(effect_term)
    if (nrow(d) == 0) {
      plot.new()
      text(
        0.5,
        0.5,
        "Run at least one community WP3 design to show this effect."
      )
      return(invisible(NULL))
    }

    p <- ggplot(
      d,
      aes(
        x = design_variant,
        y = .data[[y_var]],
        fill = design_variant
      )
    )

    if (!is.null(hline)) {
      p <- p +
        geom_hline(
          yintercept = hline,
          linetype = "dashed",
          color = "grey40",
          linewidth = 0.7
        )
      if (!is.null(hline_label)) {
        p <- p +
          annotate(
            "text",
            x = d$design_variant[1],
            y = hline + 0.02,
            label = hline_label,
            hjust = 0,
            size = 3.2,
            color = "grey40"
          )
      }
    }

    p <- p +
      geom_col(alpha = 0.86, width = 0.62) +
      scale_fill_manual(
        values = c("Post-only" = "#1f77b4", "Panel" = "#d62728"),
        name = "WP3 community design"
      ) +
      labs(
        x = NULL,
        y = y_label,
        title = paste(
          title_prefix,
          WP3_COMMUNITY_EFFECT_LABELS[[effect_term]]
        )
      ) +
      base_theme()

    if (!is.null(y_limits) || percent_y) {
      p <- p +
        scale_y_continuous(
          labels = if (percent_y) percent_format(accuracy = 1) else waiver()
        )
    }
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }

    p
  }

  register_custom_wp3_community_effect_plots <- function(
    prefix,
    y_var,
    y_label,
    title_prefix,
    hline = NULL,
    hline_label = NULL,
    y_limits = NULL,
    percent_y = FALSE
  ) {
    for (effect_term in WP3_COMMUNITY_EFFECT_ORDER) {
      local({
        effect_local <- effect_term
        output[[effect_output_id(prefix, effect_local)]] <- renderPlot({
          custom_wp3_community_effect_plot(
            effect_term = effect_local,
            y_var = y_var,
            y_label = y_label,
            title_prefix = title_prefix,
            hline = hline,
            hline_label = hline_label,
            y_limits = y_limits,
            percent_y = percent_y
          )
        })
      })
    }
  }

  register_custom_wp3_community_effect_plots(
    "custom_community_power",
    "power",
    "Power",
    "Power:",
    hline = 0.80,
    hline_label = "80% target",
    y_limits = c(0, 1),
    percent_y = TRUE
  )
  register_custom_wp3_community_effect_plots(
    "custom_community_coverage",
    "coverage",
    "95% CI coverage",
    "Coverage:",
    hline = 0.95,
    hline_label = "95%",
    y_limits = c(0.80, 1),
    percent_y = TRUE
  )
  register_custom_wp3_community_effect_plots(
    "custom_community_bias",
    "bias",
    "Bias (estimate - estimand)",
    "Bias:",
    hline = 0
  )
  register_custom_wp3_community_effect_plots(
    "custom_community_rmse",
    "rmse",
    "RMSE",
    "RMSE:"
  )

  output$custom_memory_status <- renderUI({
    memory <- custom_memory()
    if (nrow(memory) == 0) {
      return(div(
        class = "alert alert-info",
        icon("info-circle"),
        " Completed simulations will be appended here and saved as CSV files in ",
        tags$code(custom_memory_dir),
        "."
      ))
    }

    div(
      class = "alert alert-secondary",
      icon("table"),
      " Stored ",
      length(unique(memory$simulation_id)),
      " simulation",
      if (length(unique(memory$simulation_id)) == 1) "" else "s",
      " across ",
      nrow(memory),
      " wave/design row",
      if (nrow(memory) == 1) "" else "s",
      ". CSV files are saved in ",
      tags$code(custom_memory_dir),
      "."
    )
  })

  output$custom_memory_table <- renderDT({
    memory <- custom_memory()

    display <- memory |>
      rename(
        "Simulation ID" = simulation_id,
        "Saved at" = saved_at,
        "CSV file" = csv_file,
        "Waves" = waves,
        "Design type" = design_type,
        "N" = N,
        "tau" = tau,
        "J" = J,
        "Flooded community share" = p_comm_flooded,
        "Individual flooding rate" = p_indiv_given_comm,
        "Community flooding effect" = tau_community,
        "Individual flooding effect" = tau_individual,
        "rho" = rho_y,
        "Attrition" = attrition_rate,
        "Survey retention penalty" = differential_attrition,
        "Survey assignment prob" = treat_prob,
        "Flood retention boost" = flood_response_boost,
        "Sims" = sims,
        "Power: survey" = power_Z,
        "Power: community flooding" = power_community_flooded,
        "Power: individual flooding" = power_individual_flooded,
        "Power: total direct damage" = power_total_direct_effect_latent,
        "Power: survey x community" = `power_Z:community_flooded`,
        "Power: survey x individual" = `power_Z:individual_flooded`
      )

    datatable(
      display,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list("csv", "excel"),
        scrollX = TRUE,
        order = list(list(0, "desc"), list(3, "asc")),
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      )
    ) |>
      formatRound(
        c(
          "tau",
          "J",
          "Flooded community share",
          "Individual flooding rate",
          "Community flooding effect",
          "Individual flooding effect",
          "rho",
          "Attrition",
          "Survey retention penalty",
          "Survey assignment prob",
          "Flood retention boost",
          "Power: survey",
          "Power: community flooding",
          "Power: individual flooding",
          "Power: total direct damage",
          "Power: survey x community",
          "Power: survey x individual"
        ),
        digits = 3
      )
  })

  output$custom_table <- renderDT({
    r <- custom_results()
    req(!is.null(r))
    r |>
      select(
        design_label_display,
        estimator_label,
        term,
        inquiry,
        outcome_scale,
        N,
        tau,
        J,
        p_comm_flooded,
        p_indiv_given_comm,
        tau_community,
        tau_individual,
        rho_y,
        attrition_rate,
        differential_attrition,
        treat_prob,
        flood_response_boost,
        power,
        bias,
        rmse,
        coverage,
        mean_se,
        type_s_error,
        mean_estimate,
        mean_estimand
      ) |>
      rename(
        "Design" = design_label_display,
        "Estimator" = estimator_label,
        "Effect" = term,
        "Inquiry" = inquiry,
        "Scale" = outcome_scale,
        "N" = N,
        "tau (SD)" = tau,
        "J" = J,
        "Flooded community share" = p_comm_flooded,
        "Individual flooding rate" = p_indiv_given_comm,
        "Community flooding effect" = tau_community,
        "Individual flooding effect" = tau_individual,
        "rho" = rho_y,
        "Attrition" = attrition_rate,
        "Survey retention penalty" = differential_attrition,
        "Survey assignment prob" = treat_prob,
        "Flood retention boost" = flood_response_boost,
        "Power" = power,
        "Bias" = bias,
        "RMSE" = rmse,
        "Coverage" = coverage,
        "Mean SE" = mean_se,
        "Type S" = type_s_error,
        "Mean estimate" = mean_estimate,
        "Mean estimand" = mean_estimand
      ) |>
      datatable(
        rownames = FALSE,
        options = list(pageLength = 20, scrollX = TRUE, dom = "t")
      ) |>
      formatRound(
        c(
          "Power",
          "Bias",
          "RMSE",
          "Coverage",
          "Mean SE",
          "Type S",
          "Mean estimate",
          "Mean estimand"
        ),
        digits = 3
      )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Launch
# ══════════════════════════════════════════════════════════════════════════════

shinyApp(ui = ui, server = server)
