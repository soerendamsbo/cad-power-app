# design_builder.R
# Shared design components for the WP2/WP3 power analysis app.
# Sourced by precompute.R and app.R.
#
# -- HOW TO EXTEND ------------------------------------------------------------
# To add a new estimator or design type:
#   1. Add estimator declaration(s) in build_wp2_design() below.
#   2. Add a new `type` branch in the assembly block at the bottom of that
#      function.
#   3. Register the new type in DESIGN_REGISTRY.
#   4. Add a display label in ESTIMATOR_LABELS.
# No other files need to change for the new design to appear in the app.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DeclareDesign)
  library(fabricatr)
  library(randomizr)
  library(estimatr)
  library(dplyr)
})

# -- Helpers ------------------------------------------------------------------

# Map a continuous latent variable to a 1-5 Likert-style scale.
# Cut-points are chosen so that a standard-normal latent variable yields a
# roughly bell-shaped distribution across the five categories.
latent_to_likert <- function(x) {
  as.integer(cut(
    x,
    breaks = c(-Inf, -1.2, -0.35, 0.35, 1.2, Inf),
    labels = FALSE
  ))
}

# Diagnosands computed for every design cell.
# To add a diagnosand, append a line here. It will appear automatically in
# all plots and tables throughout the app.
diagnosands_power <- declare_diagnosands(
  power = mean(p.value < 0.05, na.rm = TRUE),
  bias = mean(estimate - estimand, na.rm = TRUE),
  rmse = sqrt(mean((estimate - estimand)^2, na.rm = TRUE)),
  coverage = mean(
    conf.low <= estimand & conf.high >= estimand,
    na.rm = TRUE
  ),
  mean_estimate = mean(estimate, na.rm = TRUE),
  mean_estimand = mean(estimand, na.rm = TRUE),
  mean_se = mean(std.error, na.rm = TRUE),
  type_s_error = mean(
    p.value < 0.05 & sign(estimate) != sign(estimand),
    na.rm = TRUE
  )
)

# Estimate the full set of WP3 effects from a linear interaction model. For
# linear models these are the same contrasts targeted by marginal effects:
# subgroup effects are coefficient combinations and averages weight the
# interaction by the observed distribution of the other treatment.
wp3_effects_estimator <- function(data, outcome, adjust_pre = FALSE, term = NULL) {
  rhs <- if (adjust_pre) {
    "Z * flooded + Y_pre_latent"
  } else {
    "Z * flooded"
  }

  fit <- lm(stats::as.formula(paste(outcome, "~", rhs)), data = data)
  coefs <- stats::coef(fit)
  vc <- stats::vcov(fit)
  df <- stats::df.residual(fit)
  crit <- stats::qt(0.975, df = df)

  p_flood <- mean(data$flooded, na.rm = TRUE)
  p_z <- mean(data$Z, na.rm = TRUE)

  combos <- list(
    survey_avg = c(Z = 1, `Z:flooded` = p_flood),
    survey_when_not_flooded = c(Z = 1),
    survey_when_flooded = c(Z = 1, `Z:flooded` = 1),
    survey_flood_contrast = c(`Z:flooded` = 1),
    flood_avg = c(flooded = 1, `Z:flooded` = p_z),
    flood_when_survey_control = c(flooded = 1),
    flood_when_survey_treated = c(flooded = 1, `Z:flooded` = 1),
    flood_survey_contrast = c(`Z:flooded` = 1)
  )

  estimate_combo <- function(term, weights) {
    w <- setNames(rep(0, length(coefs)), names(coefs))
    w[names(weights)] <- weights

    estimate <- sum(w * coefs)
    std.error <- sqrt(as.numeric(t(w) %*% vc %*% w))
    statistic <- estimate / std.error
    p.value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)

    data.frame(
      term = term,
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      conf.low = estimate - crit * std.error,
      conf.high = estimate + crit * std.error
    )
  }

  dplyr::bind_rows(Map(estimate_combo, names(combos), combos))
}

# Estimate WP3 community/property flooding effects with community-clustered CR2
# standard errors. The total direct-damage effect is the linear combination of
# community_flooded + individual_flooded.
wp3_community_effects_estimator <- function(
  data,
  outcome,
  adjust_pre = FALSE,
  term = NULL
) {
  out_terms <- c(
    "Z",
    "community_flooded",
    "individual_flooded",
    "total_direct_effect_latent",
    "Z:community_flooded",
    "Z:individual_flooded"
  )

  empty_result <- function() {
    data.frame(
      term = out_terms,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
  }

  rhs <- if (adjust_pre) {
    "Z * community_flooded + Z * individual_flooded + Y_pre_latent"
  } else {
    "Z * community_flooded + Z * individual_flooded"
  }

  fit <- tryCatch(
    estimatr::lm_robust(
      stats::as.formula(paste(outcome, "~", rhs)),
      data = data,
      clusters = community_id,
      se_type = "CR2"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(empty_result())
  }

  coefs <- stats::coef(fit)
  vc <- stats::vcov(fit)
  conf <- stats::confint(fit)

  direct_term <- function(term_name) {
    if (!term_name %in% names(coefs)) {
      return(data.frame(
        term = term_name,
        estimate = NA_real_,
        std.error = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      ))
    }

    data.frame(
      term = term_name,
      estimate = unname(coefs[[term_name]]),
      std.error = unname(fit$std.error[[term_name]]),
      statistic = unname(fit$statistic[[term_name]]),
      p.value = unname(fit$p.value[[term_name]]),
      conf.low = unname(conf[term_name, 1]),
      conf.high = unname(conf[term_name, 2])
    )
  }

  total_direct <- function() {
    needed <- c("community_flooded", "individual_flooded")
    if (!all(needed %in% names(coefs))) {
      return(data.frame(
        term = "total_direct_effect_latent",
        estimate = NA_real_,
        std.error = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      ))
    }

    w <- setNames(rep(0, length(coefs)), names(coefs))
    w[needed] <- 1
    estimate <- sum(w * coefs)
    std.error <- sqrt(as.numeric(t(w) %*% vc %*% w))
    statistic <- estimate / std.error
    df <- min(fit$df[needed], na.rm = TRUE)
    p.value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
    crit <- stats::qt(0.975, df = df)

    data.frame(
      term = "total_direct_effect_latent",
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      conf.low = estimate - crit * std.error,
      conf.high = estimate + crit * std.error
    )
  }

  dplyr::bind_rows(
    direct_term("Z"),
    direct_term("community_flooded"),
    direct_term("individual_flooded"),
    total_direct(),
    direct_term("Z:community_flooded"),
    direct_term("Z:individual_flooded")
  )
}

# -- Design builder ------------------------------------------------------------

#' Build a WP2/WP3 survey design.
#'
#' @param N                      Wave-1 sample size.
#' @param tau                    Average survey-embedded treatment effect in
#'                               latent SD units.
#' @param tau_flood              Average flood exposure effect in latent SD
#'                               units (typically negative for trust outcomes).
#' @param tau_interaction        Interaction effect between survey experiment and
#'                               flood exposure in latent SD units.
#' @param rho_y                  Test-retest (wave-1 to wave-2) correlation of
#'                               the outcome under no treatment. Bounded [-1,1].
#' @param attrition_rate         Baseline probability of non-response at wave 2.
#'                               Panel designs only; ignored by post-only.
#' @param differential_attrition Extra dropout log-odds for survey experiment-
#'                               treated units.
#'                               0 = no differential attrition.
#' @param flood_exposure_rate    Proportion naturally exposed to flooding.
#' @param treat_prob             Probability of survey-embedded assignment Z = 1.
#' @param flood_response_boost   Extra retention log-odds at wave 2 for flood-
#'                               exposed respondents (positive lowers attrition;
#'                               negative raises attrition).
#' @param type                   One of:
#'                               wp2_post_only, wp2_panel,
#'                               wp3_post_only, wp3_panel.
#'                               Backward-compatible aliases: post_only, panel.
#' @return A DeclareDesign object.
build_wp2_design <- function(
  N,
  tau,
  tau_flood = -0.10,
  tau_interaction = 0.04,
  rho_y,
  attrition_rate = 0.30,
  differential_attrition = 0.00,
  flood_exposure_rate = 0.10,
  treat_prob = 0.50,
  flood_response_boost = 0.00,
  type = c(
    "wp2_post_only",
    "wp2_panel",
    "wp3_post_only",
    "wp3_panel"
  )
) {
  # Backward compatibility for existing calls.
  if (length(type) == 1 && type == "post_only") {
    type <- "wp2_post_only"
  }
  if (length(type) == 1 && type == "panel") {
    type <- "wp2_panel"
  }

  type <- match.arg(type)

  # -- Model: population and potential outcomes -------------------------------
  model <- declare_model(
    N = N,
    # Demographic covariates (Danish adult population approximation)
    gender = factor(sample(
      c("Woman", "Man"),
      N,
      replace = TRUE,
      prob = c(0.51, 0.49)
    )),
    age_group = factor(sample(
      c("18-34", "35-54", "55+"),
      N,
      replace = TRUE,
      prob = c(0.28, 0.35, 0.37)
    )),
    education = factor(sample(
      c("Low", "Medium", "High"),
      N,
      replace = TRUE,
      prob = c(0.25, 0.45, 0.30)
    )),
    ideology3 = factor(sample(
      c("Left", "Center", "Right"),
      N,
      replace = TRUE,
      prob = c(0.34, 0.33, 0.33)
    )),
    # Climate concern is correlated with ideology
    climate_high = rbinom(
      N,
      1,
      ifelse(
        ideology3 == "Left",
        0.72,
        ifelse(ideology3 == "Center", 0.55, 0.38)
      )
    ),
    # Natural experiment exposure indicator (WP3)
    flood_exposed = rbinom(N, 1, flood_exposure_rate),
    # Wave-1 latent outcome: linear combination of demographics + noise
    baseline_mu = 0.20 *
      (gender == "Woman") +
      0.15 * (age_group == "55+") +
      0.10 * (education == "High") +
      0.35 * (ideology3 == "Left") -
      0.25 * (ideology3 == "Right") +
      0.45 * climate_high,
    Y_pre_raw = baseline_mu + rnorm(N),
    Y_pre_latent = as.numeric(scale(Y_pre_raw)),
    # Wave-2 control potential outcome: rho_y pins the test-retest correlation.
    # Y0_latent = rho * Y_pre + sqrt(1 - rho^2) * shock preserves unit variance.
    post_shock = rnorm(N),
    Y0_latent = rho_y * Y_pre_latent + sqrt(1 - rho_y^2) * post_shock,
    # Survey-treatment heterogeneity around tau.
    tau_z_i = tau +
      0.03 * climate_high -
      0.02 * (ideology3 == "Right"),
    # Potential outcomes for the 2x2 design: Z in {0,1}, F in {0,1}
    Y_latent_Z_0_F_0 = Y0_latent,
    Y_latent_Z_1_F_0 = Y0_latent + tau_z_i,
    Y_latent_Z_0_F_1 = Y0_latent + tau_flood,
    Y_latent_Z_1_F_1 = Y0_latent + tau_z_i + tau_flood + tau_interaction,
    # Coarsened Likert outcomes
    Y_1_5_Z_0_F_0 = latent_to_likert(Y_latent_Z_0_F_0),
    Y_1_5_Z_1_F_0 = latent_to_likert(Y_latent_Z_1_F_0),
    Y_1_5_Z_0_F_1 = latent_to_likert(Y_latent_Z_0_F_1),
    Y_1_5_Z_1_F_1 = latent_to_likert(Y_latent_Z_1_F_1)
  )

  # -- Inquiries: estimands ---------------------------------------------------
  # ## ADD NEW ESTIMANDS HERE (step 1 of 2 - also add estimator below)
  inquiry_wp2 <- declare_inquiry(
    ATE_latent = mean(Y_latent_Z_1_F_0 - Y_latent_Z_0_F_0),
    ATE_1_5 = mean(Y_1_5_Z_1_F_0 - Y_1_5_Z_0_F_0)
  )

  inquiry_wp3 <- declare_inquiry(
    survey_avg_latent = mean(
      flood_exposed * (Y_latent_Z_1_F_1 - Y_latent_Z_0_F_1) +
        (1 - flood_exposed) * (Y_latent_Z_1_F_0 - Y_latent_Z_0_F_0)
    ),
    survey_when_not_flooded_latent = mean(
      Y_latent_Z_1_F_0 - Y_latent_Z_0_F_0
    ),
    survey_when_flooded_latent = mean(
      Y_latent_Z_1_F_1 - Y_latent_Z_0_F_1
    ),
    survey_flood_contrast_latent = mean(
      (Y_latent_Z_1_F_1 - Y_latent_Z_0_F_1) -
        (Y_latent_Z_1_F_0 - Y_latent_Z_0_F_0)
    ),
    flood_avg_latent = mean(
      Z * (Y_latent_Z_1_F_1 - Y_latent_Z_1_F_0) +
        (1 - Z) * (Y_latent_Z_0_F_1 - Y_latent_Z_0_F_0)
    ),
    flood_when_survey_control_latent = mean(
      Y_latent_Z_0_F_1 - Y_latent_Z_0_F_0
    ),
    flood_when_survey_treated_latent = mean(
      Y_latent_Z_1_F_1 - Y_latent_Z_1_F_0
    ),
    flood_survey_contrast_latent = mean(
      (Y_latent_Z_1_F_1 - Y_latent_Z_1_F_0) -
        (Y_latent_Z_0_F_1 - Y_latent_Z_0_F_0)
    ),
    survey_avg_1_5 = mean(
      flood_exposed * (Y_1_5_Z_1_F_1 - Y_1_5_Z_0_F_1) +
        (1 - flood_exposed) * (Y_1_5_Z_1_F_0 - Y_1_5_Z_0_F_0)
    ),
    survey_when_not_flooded_1_5 = mean(Y_1_5_Z_1_F_0 - Y_1_5_Z_0_F_0),
    survey_when_flooded_1_5 = mean(Y_1_5_Z_1_F_1 - Y_1_5_Z_0_F_1),
    survey_flood_contrast_1_5 = mean(
      (Y_1_5_Z_1_F_1 - Y_1_5_Z_0_F_1) -
        (Y_1_5_Z_1_F_0 - Y_1_5_Z_0_F_0)
    ),
    flood_avg_1_5 = mean(
      Z * (Y_1_5_Z_1_F_1 - Y_1_5_Z_1_F_0) +
        (1 - Z) * (Y_1_5_Z_0_F_1 - Y_1_5_Z_0_F_0)
    ),
    flood_when_survey_control_1_5 = mean(Y_1_5_Z_0_F_1 - Y_1_5_Z_0_F_0),
    flood_when_survey_treated_1_5 = mean(Y_1_5_Z_1_F_1 - Y_1_5_Z_1_F_0),
    flood_survey_contrast_1_5 = mean(
      (Y_1_5_Z_1_F_1 - Y_1_5_Z_1_F_0) -
        (Y_1_5_Z_0_F_1 - Y_1_5_Z_0_F_0)
    )
  )

  # -- Data strategy ----------------------------------------------------------
  assignment <- declare_assignment(
    Z = randomizr::complete_ra(N = N, prob = treat_prob)
  )

  measurement <- declare_measurement(
    flooded = flood_exposed,
    Y_latent = ifelse(
      Z == 1 & flooded == 1,
      Y_latent_Z_1_F_1,
      ifelse(
        Z == 1 & flooded == 0,
        Y_latent_Z_1_F_0,
        ifelse(
          Z == 0 & flooded == 1,
          Y_latent_Z_0_F_1,
          Y_latent_Z_0_F_0
        )
      )
    ),
    Y_1_5 = ifelse(
      Z == 1 & flooded == 1,
      Y_1_5_Z_1_F_1,
      ifelse(
        Z == 1 & flooded == 0,
        Y_1_5_Z_1_F_0,
        ifelse(
          Z == 0 & flooded == 1,
          Y_1_5_Z_0_F_1,
          Y_1_5_Z_0_F_0
        )
      )
    )
  )

  # Attrition block (panel designs only): response probability at wave 2 is a
  # logistic function of treatment, flood exposure, climate concern, and
  # education. Positive flood_response_boost means flood exposure increases
  # follow-up participation; negative values lower follow-up participation.
  attrition <- declare_measurement(
    responds_post = rbinom(
      N,
      1,
      plogis(
        qlogis(1 - attrition_rate) -
          differential_attrition * Z +
          flood_response_boost * flooded +
          0.15 * climate_high -
          0.10 * (education == "Low")
      )
    ),
    filter = responds_post == 1
  )

  # -- Answer strategy: estimators --------------------------------------------
  wp2_post_only_estimators <-
    declare_estimator(
      Y_latent ~ Z,
      .method = lm,
      inquiry = "ATE_latent",
      term = "Z",
      label = "wp2_post_only_latent"
    ) +
    declare_estimator(
      Y_1_5 ~ Z,
      .method = lm,
      inquiry = "ATE_1_5",
      term = "Z",
      label = "wp2_post_only_1_5"
    )

  wp2_panel_estimators <-
    declare_estimator(
      Y_latent ~ Z + Y_pre_latent,
      .method = lm,
      inquiry = "ATE_latent",
      term = "Z",
      label = "wp2_panel_adjusted_latent"
    ) +
    declare_estimator(
      Y_1_5 ~ Z + Y_pre_latent,
      .method = lm,
      inquiry = "ATE_1_5",
      term = "Z",
      label = "wp2_panel_adjusted_1_5"
    )

  wp3_latent_inquiries <- c(
    "survey_avg_latent",
    "survey_when_not_flooded_latent",
    "survey_when_flooded_latent",
    "survey_flood_contrast_latent",
    "flood_avg_latent",
    "flood_when_survey_control_latent",
    "flood_when_survey_treated_latent",
    "flood_survey_contrast_latent"
  )

  wp3_1_5_inquiries <- c(
    "survey_avg_1_5",
    "survey_when_not_flooded_1_5",
    "survey_when_flooded_1_5",
    "survey_flood_contrast_1_5",
    "flood_avg_1_5",
    "flood_when_survey_control_1_5",
    "flood_when_survey_treated_1_5",
    "flood_survey_contrast_1_5"
  )

  wp3_effect_terms <- c(
    "survey_avg",
    "survey_when_not_flooded",
    "survey_when_flooded",
    "survey_flood_contrast",
    "flood_avg",
    "flood_when_survey_control",
    "flood_when_survey_treated",
    "flood_survey_contrast"
  )

  wp3_post_only_estimators <-
    declare_estimator(
      handler = label_estimator(wp3_effects_estimator),
      outcome = "Y_latent",
      adjust_pre = FALSE,
      inquiry = wp3_latent_inquiries,
      term = wp3_effect_terms,
      label = "wp3_post_latent"
    ) +
    declare_estimator(
      handler = label_estimator(wp3_effects_estimator),
      outcome = "Y_1_5",
      adjust_pre = FALSE,
      inquiry = wp3_1_5_inquiries,
      term = wp3_effect_terms,
      label = "wp3_post_1_5"
    )

  wp3_panel_estimators <-
    declare_estimator(
      handler = label_estimator(wp3_effects_estimator),
      outcome = "Y_latent",
      adjust_pre = TRUE,
      inquiry = wp3_latent_inquiries,
      term = wp3_effect_terms,
      label = "wp3_panel_latent"
    ) +
    declare_estimator(
      handler = label_estimator(wp3_effects_estimator),
      outcome = "Y_1_5",
      adjust_pre = TRUE,
      inquiry = wp3_1_5_inquiries,
      term = wp3_effect_terms,
      label = "wp3_panel_1_5"
    )

  # ## ADD NEW ESTIMATORS HERE (step 2 of 2)

  # -- Assemble design --------------------------------------------------------
  # ## ADD NEW DESIGN TYPES HERE: add an else-if branch and register below.
  if (type == "wp2_post_only") {
    model + inquiry_wp2 + assignment + measurement + wp2_post_only_estimators
  } else if (type == "wp2_panel") {
    model + inquiry_wp2 + assignment + measurement + attrition + wp2_panel_estimators
  } else if (type == "wp3_post_only") {
    model + assignment + inquiry_wp3 + measurement + wp3_post_only_estimators
  } else if (type == "wp3_panel") {
    model + assignment + inquiry_wp3 + measurement + attrition + wp3_panel_estimators
  }
}

#' Build a WP3 community/property flooding design.
#'
#' @param N                      Total sample size.
#' @param J                      Number of communities.
#' @param p_comm_flooded         Proportion of communities that flood.
#' @param p_indiv_given_comm     Proportion directly flooded within flooded
#'                               communities.
#' @param tau                    Average survey-embedded treatment effect in
#'                               latent SD units.
#' @param tau_community          Community-level flooding effect in latent SD
#'                               units.
#' @param tau_individual         Incremental individual/property flooding effect
#'                               within flooded communities.
#' @param rho_y                  Wave-1 to wave-2 correlation.
#' @param attrition_rate         Baseline probability of non-response at wave 2.
#' @param differential_attrition Extra dropout log-odds for survey-treated units.
#' @param treat_prob             Survey experiment assignment probability.
#' @param flood_response_boost   Extra retention log-odds for community-flooded
#'                               respondents.
#' @param type                   wp3_community_post_only or wp3_community_panel.
#' @return A DeclareDesign object.
build_wp3_community_design <- function(
  N,
  J = 200,
  p_comm_flooded = 0.15,
  p_indiv_given_comm = 0.40,
  tau = 0.10,
  tau_community = -0.10,
  tau_individual = -0.05,
  rho_y,
  attrition_rate = 0.30,
  differential_attrition = 0.00,
  treat_prob = 0.50,
  flood_response_boost = 0.00,
  type = c(
    "wp3_community_post_only",
    "wp3_community_panel"
  )
) {
  type <- match.arg(type)

  model <- declare_model(
    N = N,
    community_id = sample(seq_len(J), N, replace = TRUE),
    community_flooded = {
      community_flooded_j <- stats::rbinom(J, 1, p_comm_flooded)
      community_flooded_j[community_id]
    },
    individual_flooded = community_flooded *
      stats::rbinom(N, 1, p_indiv_given_comm),
    gender = factor(sample(
      c("Woman", "Man"),
      N,
      replace = TRUE,
      prob = c(0.51, 0.49)
    )),
    age_group = factor(sample(
      c("18-34", "35-54", "55+"),
      N,
      replace = TRUE,
      prob = c(0.28, 0.35, 0.37)
    )),
    education = factor(sample(
      c("Low", "Medium", "High"),
      N,
      replace = TRUE,
      prob = c(0.25, 0.45, 0.30)
    )),
    ideology3 = factor(sample(
      c("Left", "Center", "Right"),
      N,
      replace = TRUE,
      prob = c(0.34, 0.33, 0.33)
    )),
    climate_high = rbinom(
      N,
      1,
      ifelse(
        ideology3 == "Left",
        0.72,
        ifelse(ideology3 == "Center", 0.55, 0.38)
      )
    ),
    baseline_mu = 0.20 *
      (gender == "Woman") +
      0.15 * (age_group == "55+") +
      0.10 * (education == "High") +
      0.35 * (ideology3 == "Left") -
      0.25 * (ideology3 == "Right") +
      0.45 * climate_high,
    Y_pre_raw = baseline_mu + rnorm(N),
    Y_pre_latent = as.numeric(scale(Y_pre_raw)),
    post_shock = rnorm(N),
    Y0_latent = rho_y * Y_pre_latent + sqrt(1 - rho_y^2) * post_shock,
    tau_z_i = tau +
      0.03 * climate_high -
      0.02 * (ideology3 == "Right"),
    Y_latent_Z0_C0 = Y0_latent,
    Y_latent_Z1_C0 = Y0_latent + tau_z_i,
    Y_latent_Z0_C1_I0 = Y0_latent + tau_community,
    Y_latent_Z1_C1_I0 = Y0_latent + tau_z_i + tau_community,
    Y_latent_Z0_C1_I1 = Y0_latent + tau_community + tau_individual,
    Y_latent_Z1_C1_I1 = Y0_latent + tau_z_i + tau_community + tau_individual
  )

  inquiry <- declare_inquiry(
    community_effect_latent = mean(Y_latent_Z0_C1_I0 - Y_latent_Z0_C0),
    individual_effect_latent = mean(Y_latent_Z0_C1_I1 - Y_latent_Z0_C1_I0),
    survey_effect_latent = mean(Y_latent_Z1_C0 - Y_latent_Z0_C0),
    total_direct_effect_latent = mean(Y_latent_Z0_C1_I1 - Y_latent_Z0_C0)
  )

  assignment <- declare_assignment(
    Z = randomizr::complete_ra(N = N, prob = treat_prob)
  )

  measurement <- declare_measurement(
    Y_latent = ifelse(
      community_flooded == 0 & Z == 0,
      Y_latent_Z0_C0,
      ifelse(
        community_flooded == 0 & Z == 1,
        Y_latent_Z1_C0,
        ifelse(
          individual_flooded == 0 & Z == 0,
          Y_latent_Z0_C1_I0,
          ifelse(
            individual_flooded == 0 & Z == 1,
            Y_latent_Z1_C1_I0,
            ifelse(
              Z == 0,
              Y_latent_Z0_C1_I1,
              Y_latent_Z1_C1_I1
            )
          )
        )
      )
    )
  )

  attrition <- declare_measurement(
    responds_post = rbinom(
      N,
      1,
      plogis(
        qlogis(1 - attrition_rate) -
          differential_attrition * Z +
          flood_response_boost * community_flooded +
          0.15 * climate_high -
          0.10 * (education == "Low")
      )
    ),
    filter = responds_post == 1
  )

  community_terms <- c(
    "Z",
    "community_flooded",
    "individual_flooded",
    "total_direct_effect_latent",
    "Z:community_flooded",
    "Z:individual_flooded"
  )

  community_inquiries <- c(
    "survey_effect_latent",
    "community_effect_latent",
    "individual_effect_latent",
    "total_direct_effect_latent",
    NA_character_,
    NA_character_
  )

  post_only_estimators <- declare_estimator(
    handler = label_estimator(wp3_community_effects_estimator),
    outcome = "Y_latent",
    adjust_pre = FALSE,
    inquiry = community_inquiries,
    term = community_terms,
    label = "wp3_community_post_latent"
  )

  panel_estimators <- declare_estimator(
    handler = label_estimator(wp3_community_effects_estimator),
    outcome = "Y_latent",
    adjust_pre = TRUE,
    inquiry = community_inquiries,
    term = community_terms,
    label = "wp3_community_panel_latent"
  )

  if (type == "wp3_community_post_only") {
    model + assignment + inquiry + measurement + post_only_estimators
  } else if (type == "wp3_community_panel") {
    model + assignment + inquiry + measurement + attrition + panel_estimators
  }
}

# -- Registries ---------------------------------------------------------------

# Register every design type that should appear in the app.
# uses_attrition = TRUE causes attrition controls to be active for that type.
DESIGN_REGISTRY <- list(
  list(
    type = "wp2_post_only",
    label = "WP2 - Post-only (no baseline adjustment)",
    uses_attrition = FALSE
  ),
  list(
    type = "wp2_panel",
    label = "WP2 - Panel (baseline-adjusted, with attrition)",
    uses_attrition = TRUE
  ),
  list(
    type = "wp3_post_only",
    label = "WP3 - Post-only (survey + flood + interaction)",
    uses_attrition = FALSE
  ),
  list(
    type = "wp3_panel",
    label = "WP3 - Panel (survey + flood + interaction)",
    uses_attrition = TRUE
  ),
  list(
    type = "wp3_community_post_only",
    label = "WP3 - Community + property flooding (post-only)",
    uses_attrition = FALSE
  ),
  list(
    type = "wp3_community_panel",
    label = "WP3 - Community + property flooding (panel)",
    uses_attrition = TRUE
  )
)

# Human-readable names for estimators in plots and tables.
ESTIMATOR_LABELS <- c(
  wp2_post_only_latent = "WP2 Post-only - Survey experiment (latent SD)",
  wp2_post_only_1_5 = "WP2 Post-only - Survey experiment (Likert 1-5)",
  wp2_panel_adjusted_latent = "WP2 Panel - Survey experiment (latent SD)",
  wp2_panel_adjusted_1_5 = "WP2 Panel - Survey experiment (Likert 1-5)",
  wp3_post_latent = "WP3 Post-only - Effects (latent SD)",
  wp3_post_1_5 = "WP3 Post-only - Effects (Likert 1-5)",
  wp3_panel_latent = "WP3 Panel - Effects (latent SD)",
  wp3_panel_1_5 = "WP3 Panel - Effects (Likert 1-5)",
  wp3_community_post_latent = "WP3 Community - Post-only (latent SD)",
  wp3_community_panel_latent = "WP3 Community - Panel (latent SD)"
)

WP3_COMMUNITY_EFFECT_ORDER <- c(
  "Z",
  "community_flooded",
  "individual_flooded",
  "total_direct_effect_latent",
  "Z:community_flooded",
  "Z:individual_flooded"
)

WP3_COMMUNITY_EFFECT_LABELS <- c(
  community_flooded = "Community flooding",
  individual_flooded = "Individual flooding",
  Z = "Survey experiment",
  total_direct_effect_latent = "Total direct damage",
  `Z:community_flooded` = "Survey x community flooding",
  `Z:individual_flooded` = "Survey x individual flooding"
)

WP3_COMMUNITY_EFFECT_X <- c(
  community_flooded = "tau_community",
  individual_flooded = "tau_individual",
  Z = "tau",
  total_direct_effect_latent = "tau_individual",
  `Z:community_flooded` = "tau",
  `Z:individual_flooded` = "tau"
)

WP3_COMMUNITY_EFFECT_X_LABELS <- c(
  tau = "Survey experiment effect tau (SD units)",
  tau_community = "Community flooding effect tau_community (SD units)",
  tau_individual = "Individual flooding effect tau_individual (SD units)"
)
