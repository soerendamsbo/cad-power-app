# Spatial Clustering, Statistical Power, and Standard Errors in Flood Exposure Research

## Context
Research design for estimating the causal effect of flood exposure in coastal Denmark. N ≈ 4,000. Two-wave panel design with pretest values. Treatment (flood exposure) is geographically clustered. Community-level vs. individual-level exposure definitions considered.

---

## Q1: How does geographic clustering of treatment affect power? What's the difference between ordinary robust SEs, cluster-robust SEs, and Conley SEs?

### The core problem: what's your effective sample size?

Start with an extreme thought experiment. Suppose all your treated individuals live in one city and all your untreated individuals live in another. You observe 10,000 people in each. Your estimated treatment effect is just the difference in mean outcomes across the two cities. But you don't really have 20,000 independent observations — you essentially have two data points: one treated cluster and one untreated cluster. Any difference you find could just be a city-level idiosyncrasy that has nothing to do with flooding. You have no power at all.

Now relax the extreme. Suppose treatment is spread across, say, 15 municipalities, and control is spread across another 15. That's much better, but your effective sample size for estimating the treatment effect is still closer to 30 than to whatever your individual N is. Adding more individuals within those municipalities helps you estimate each municipality's mean more precisely, but it does almost nothing to help you distinguish the treatment effect from municipality-level confounds or correlated shocks.

**Power for your treatment effect depends primarily on the number of independent geographic units across which treatment varies, not on the total number of individuals.** The individual-level N buys you precision within clusters, but what you need is variation *across* clusters.

### Why treatment clustering specifically matters

Spatial correlation in outcomes alone isn't catastrophic — it inflates standard errors but can be handled. What makes this setting particularly challenging is the combination of two things: outcomes are spatially correlated (nearby people have similar economic, health, and social outcomes for reasons unrelated to flooding) *and* treatment is spatially clustered (nearby people share the same treatment status). When both are true, the bias in naive standard errors is most severe.

The intuition is that spatial correlation in the errors eats into your effective sample size, but it only matters for inference about a coefficient to the extent that the *regressor of interest* also varies along the same spatial dimension. If treatment were randomly sprinkled across space at the individual level, spatial error correlation would matter less for your treatment effect estimate because the correlation in errors wouldn't align systematically with treatment status. But floods don't work that way — they blanket entire areas — so the error correlation and the treatment variation are riding on the same geographic structure.

### The three approaches

**Ordinary robust (HC) standard errors** assume observations are independent. They correct for heteroskedasticity but nothing else. In this setting, they will be far too small because they count every individual as providing independent information about the treatment effect. If you have 5,000 treated people in a flood zone that spans three neighborhoods, HC standard errors treat that as 5,000 independent pieces of evidence. In reality, once you know outcomes for a few hundred people in that area, the additional thousands are largely redundant for identifying the treatment effect. The result is dramatic over-rejection — you'll find "significant" effects far more often than your nominal rate suggests.

**Cluster-robust standard errors** allow arbitrary within-cluster correlation but assume independence across clusters. This is much better, but two things matter enormously.

First, *the level at which you cluster.* Clustering at a small geographic level (say, census tracts) allows for correlation within tracts but still assumes independence across nearby tracts. If the spatial correlation in outcomes extends beyond tract boundaries — which it almost certainly does for anything related to flooding, local labor markets, infrastructure, and so on — you're still understating uncertainty. Clustering at a coarser level (municipalities, counties) is typically more conservative and more appropriate, because you want to cluster at or above the level at which treatment is effectively "assigned." Since floods hit contiguous areas that often span multiple small administrative units, a coarser clustering level better captures the true correlation structure.

Second, *the number of clusters.* Cluster-robust standard errors rely on the number of clusters being large for their asymptotic properties to hold. With fewer than roughly 40–50 clusters, they can be badly biased (usually downward, making you overconfident again). If your flood-affected region only gives you, say, 12 treated municipalities and 15 control municipalities, you're in the zone where cluster-robust inference is shaky and you may need small-sample corrections (like the wild cluster bootstrap).

**Conley standard errors** take a different approach. Instead of putting observations into discrete bins, they allow errors to be correlated as a continuous function of distance. You specify a kernel and a bandwidth: observations within the bandwidth distance can have correlated errors, with the correlation decaying as distance increases. This is appealing because spatial correlation doesn't actually stop at administrative boundaries — a census tract boundary doesn't make the people on either side of it suddenly independent. Conley SEs respect the continuous geography.

The tradeoff is that you have to choose the bandwidth, and results can be sensitive to this choice. Too small a bandwidth and you're back to underestimating standard errors; too large and you may be absorbing correlation that isn't really there (though this tends to be conservative, which is the safer direction). A common practice is to report results across a range of bandwidths to show robustness.

### Relative magnitudes

Roughly, you should expect the following ordering of standard errors:

HC robust < cluster-robust (fine geography) < cluster-robust (coarse geography) ≈ Conley (with appropriate bandwidth)

The jump from HC to cluster-robust can easily be a factor of 2–5x in your standard errors, which translates to a 4–25x difference in the sample size you'd need to achieve the same power.

### Recommendation

Use cluster-robust standard errors at the coarsest geographic level that still leaves you with enough clusters for reliable inference, and report Conley standard errors as a complement.

Then also report Conley SEs across a few bandwidths (say, 10km, 25km, 50km). If the Conley SEs at a reasonable bandwidth are similar to your cluster-robust SEs, that's reassuring. If they're substantially larger, that tells you spatial correlation extends beyond your cluster boundaries.

Never report only HC standard errors in your main results.

### Power, concretely

The design effect formula gives you a multiplier: DEFF = 1 + (m̄ − 1) × ρ, where m̄ is the average cluster size and ρ is the intra-cluster correlation of the outcome. Your effective sample size is N/DEFF. If you have 10,000 people in 25 clusters of 400, and the ICC is even 0.05, your DEFF is about 21, giving you an effective N of roughly 480. If the ICC is 0.10, DEFF is about 41, effective N around 245.

**The most valuable thing you can do for power is to increase the number of independent geographic clusters, not the number of individuals within clusters.**

---

## Q2: What about the definition of "exposure"? Physical property exposure vs. community exposure. Context is coastal Denmark, N ≈ 4,000, 98 municipalities.

### The treatment definition changes everything about the clustering problem

With **physical property exposure**, treatment genuinely varies at a very fine spatial scale. If one house floods and the neighbor's doesn't — due to micro-elevation, drainage, where exactly the water reached — then treatment varies nearly independently across nearby units. You could even consider within-community comparisons (community fixed effects) to absorb all community-level confounders. That's a strong design with good power properties.

With **community exposure**, treatment is definitionally constant within whatever you call a "community." Every person in the community shares the same treatment status. Individual-level variation contributes nothing to identifying the treatment effect — it only helps estimate each community's mean outcome more precisely. Your effective sample size collapses to the number of communities. This is the Moulton problem in its purest form.

### What this means for N = 4,000

With 4,000 individuals and community-level treatment, power depends on:

- How many communities do you have? If ~100–200 people per community, you're looking at 20–40 communities. If ~30–50, maybe 80–130.
- Of those, how many are treated? Coastal flooding will only hit some communities.
- How much do communities vary in outcomes for reasons unrelated to flooding?
- Within-community sample size helps only insofar as it reduces estimation error for each community's mean, with diminishing returns.

### The boundary problem

The community definition simultaneously determines: who counts as treated, how many clusters you have, and how large each cluster is.

- Broadly defined communities (parish/postal code): fewer, larger communities → less power, treatment may be diluted.
- Narrowly defined communities (a few blocks): more clusters (better for power), treatment more intense, but less plausible as social units.

Danish coastal settlements have distinct physical structure — small towns, fishing villages — so there may be natural definitions that are both socially meaningful and analytically tractable.

### Standard error choice given community-level treatment

- **HC robust SEs**: essentially meaningless. Don't use them.
- **Cluster-robust SEs**: natural baseline if you can define communities. Need 40+ clusters; otherwise use wild cluster bootstrap.
- **Conley SEs**: particularly appealing because the community boundary is fuzzy. Sidesteps boundary problem for inference (but not for treatment definition).

### The deeper tension

Property exposure gives a stronger design but captures the wrong estimand if you care about the social experience of disruption. Community exposure captures the right thing at steep statistical cost.

Consider: define community exposure as primary treatment, but also measure property exposure within treated communities. This lets you estimate the total community effect and then decompose it by direct damage — the between-community comparison gives the policy-relevant quantity, while within-community variation gives mechanism and precision.

---

## Q3: What about a highly targeted sampling strategy across many high-risk communities (100+), defined by zip code? Minority treated. N = 4,000.

### 100+ communities changes the inference picture dramatically

With 100+ clusters, cluster-robust standard errors are on very solid asymptotic footing. The few-clusters problem disappears. Straightforward cluster-robust SEs at the zip code level should work well.

Zip codes as community definition are defensible for Denmark — small enough to be socially meaningful, standardized enough to not look arbitrary.

### The binding constraint: number of treated communities

With 100+ communities but only a minority treated (say 20–30), the 1/n_treated term dominates precision. The extra control communities barely help beyond a certain point. Effective power is largely determined by the treated communities.

### Within-community thinning is less painful than it looks

With ICC of 0.05 and m = 40 per community, per-community noise ≈ 0.074σ². With m = 200, it's ≈ 0.055σ². Quintupling within-community sample only buys ~25% variance reduction. At 30–40 per community, you've captured most of the available gains.

### Allocation strategy

Oversample treated communities. If 25 treated and 85 control communities: instead of ~36 per community everywhere, consider ~80 per treated community and ~22 per control community. This sharpens treated community means and partially compensates for cluster count imbalance.

### Targeted sampling helps power

High-risk coastal communities may be more similar to each other than a random sample of Danish communities, reducing σ²_between. Targeted sampling into a more homogeneous population is a feature.

---

## Q4: What if we sample only a few people per community, spreading across many communities? Use covariates and pretest values to recoup precision.

### Why thin spreading works

For fixed N, you almost always gain more power by maximizing number of clusters than within-cluster sample size.

**Comparison with N = 4,000 and 30% treated:**

- Design A: 100 communities × 40 people. Per-community noise ≈ 0.074σ² (ICC=0.05). Cluster-count term ≈ 0.048. Variance ∝ 0.0036.
- Design B: 500 communities × 8 people. Per-community noise ≈ 0.169σ². Cluster-count term ≈ 0.0095. Variance ∝ 0.0016.

**Design B has less than half the variance** despite each community being estimated much less precisely. Noisy community means are fine as long as that noise is independent across communities.

### The pretest makes thin spreading even more attractive

Controlling for Y_pre in an ANCOVA reduces residual variance by roughly (1 − ρ²). At pre-post correlation of 0.7, that's ~50% reduction. Critically, the pretest absorbs both between- and within-community variance, driving the residual ICC way down. This makes the penalty for few people per community nearly vanish.

### Practical design

With N = 4,000: aim for 400–600 zip codes with 7–10 people each. Oversample flood-likely areas to maximize treated community count (100–180 treated communities instead of 20–30).

### What you give up

- Can't do within-community analysis with 8 people per community
- Can't use community fixed effects (but ANCOVA with pretest substitutes well)
- More reliant on pretest being a good proxy
- Logistics of sampling a few people from many zip codes (but Danish register infrastructure makes this viable)

---

## Q5: Does the bottom line hold with very small clusters (1–5)? What about complete random assignment vs. cluster assignment?

### OLS with variable and very small cluster sizes

OLS is fine. With community-level treatment, OLS estimates a weighted average of community-level effects (larger communities get more weight). Valid estimand.

Cluster-robust SEs with singletons and very small clusters are not a problem when you have many clusters overall. A singleton contributes its individual residual to the sandwich estimator — valid, just noisy.

**Key limiting property:** as cluster sizes shrink toward 1, cluster-robust SEs converge toward HC robust SEs. With all singletons, they're identical. Thin spreading gradually makes the clustering correction irrelevant by design.

### Cluster assignment vs. individual assignment

**Cluster-assigned** (community exposure — everyone in flooded community is treated): treatment effect identified purely from between-community comparisons. Effective N ≈ number of communities.

**Individually assigned** (exposure varies within communities, as-if random): treatment effect identified from within-community comparisons. Community fixed effects absorb all community-level confounders. Effective N approaches number of individuals.

### The magnitude of the difference

- Individual assignment: effective N ≈ 2,500–3,500 (out of 4,000, with minor clustering adjustment)
- Cluster assignment with thin spreading (400+ communities): effective N ≈ 400–600
- Cluster assignment with fewer communities (80–100): effective N can drop below 100

**The same 4,000 people can yield effective N from under 100 to nearly 4,000** depending on treatment assignment level. That's the difference between a well-powered study and one that can only detect very large effects.

### Strategy

Define community exposure as primary treatment (between-community comparison), but use individual-level variation within treated communities for dose-response/mechanisms. First analysis = headline with cluster-appropriate inference. Second analysis = within-community identification with much more power.

---

## Q6: Can I use a standard power calculator with effective N?

### Individual-level assignment

Standard calculator works directly. Plug in N ≈ 4,000, use residual SD after pretest adjustment (multiply outcome SD by √(1 − ρ²)), done.

### Cluster-level assignment

**Cleanest approach:** treat it as a two-sample comparison where each observation is a community.

- N = number of treated communities in one arm, control communities in the other
- SD = SD of cluster means = √(σ²_between + σ²_within/m)
- Adjust for pretest: multiply by √(1 − ρ²)

The effective-N shortcut (N/DEFF with individual SD) gives approximately the same answer but is error-prone — people often double-correct.

A **cluster-RCT power calculator** is less error-prone: it asks for number of clusters per arm, individuals per cluster, ICC, and effect size directly.

### SD of cluster means vs. individual SD

With thin spreading (small m), the SD of cluster means approaches the individual SD:
- m = 1: identical
- m = 3–5: same neighborhood (~40–70% of individual SD depending on ICC)
- m = 10+: substantially smaller

---

## Q7: How to estimate cluster SD from related data in R

### Variance decomposition approach (recommended)

```r
library(lme4)

mod <- lmer(outcome ~ 1 + (1 | zipcode), data = df)
vc <- as.data.frame(VarCorr(mod))

s2_b <- vc$vcov[vc$grp == "zipcode"]
s2_w <- vc$vcov[vc$grp == "Residual"]
icc <- s2_b / (s2_b + s2_w)

# SD of cluster means for your planned m
m <- 6
sd_cluster <- sqrt(s2_b + s2_w / m)

# Adjust for pretest
rho_pre_post <- cor(df$outcome_pre, df$outcome_post)
sd_cluster_adj <- sd_cluster * sqrt(1 - rho_pre_post^2)
```

This is more useful than `sd(cluster_means)` directly because it lets you plug in different values of m and see how power changes with thinner/thicker spreading.

**Note:** If the related data covers all of Denmark, restrict to coastal zip codes before fitting the model for a more relevant ICC estimate.
