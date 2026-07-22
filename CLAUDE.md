# Sierra Gating Project — Context for Claude Code

## What this project is
Research thesis: **group-aware stochastic gating for correlated gene selection** in
single-cell transcriptomic data (scRNA-seq). Target diseases: IPF, COPD, FHP.

## The core problem
Standard feature selection (LASSO/L1, standard stochastic gating) treats every gene
as independent. Real disease biology is driven by **correlated gene modules**, not
single genes. These methods tend to pick 1-2 "representative" genes from a
correlated module and drop the rest — even if the whole module is truly relevant.

## The novel contribution
Extend stochastic gating (each gene gets a learnable "gate" trained jointly with the
classifier via gradient descent) to be **group-aware**: couple gates within known/
inferred correlated gene modules, so a whole module turns on/off together instead of
being scored gene-by-gene. Selection happens *during* training, not before/after.

## Research questions
1. Do standard L1/stochastic gating fail to recover correlated gene modules?
2. Does incorporating correlation improve feature-selection stability/reproducibility?
3. Does group-aware gating improve interpretability without sacrificing accuracy?
4. Does it generalize across diseases/simulation settings?
5. (secondary) Does multi-view integration across cell types help patient representations?

## Methodology (7 steps — see research plan for full detail)
1. Generate synthetic scRNA-seq data (Negative Binomial counts, patient/cell/gene
   hierarchy) — **done**, see `edits_syntheticdata.ipynb`.
2. Construct correlated gene modules — **done**, see "Current state" below.
3. Implement baselines: standard stochastic gates, L1, Elastic Net — **L1 done**
   (now evaluated per-module, not just signal-vs-background).
4. Build the group-aware gating framework — **not started**, the main contribution
   and the next real task.
5. Sweep simulation params (signal strength, correlation, noise, dimensionality) to
   stress-test robustness.
6. Validate on real IPF/COPD/FHP datasets.
7. Compare methods on AUROC/F1, feature recovery, stability, pathway enrichment.

## Current state of `edits_syntheticdata.ipynb`
Simulator: imports/seed → parameters (incl. module knobs) → patients_df → cells_df
(size_factor) → genes_df (signal/background split + `module_id` +
`is_disease_module`) → gamma_i (patient random effect) → group-specific means →
**module-aware Gaussian copula NB simulation** → pivot to expression_matrix →
summary stats → box plot → extended data-QC/report section (dataset overview
table, size-factor/gamma_i distributions, expression histograms, per-module
boxplots, correlation heatmaps, PCA). Followed by an L1 logistic regression
baseline test (train/test split by patient via `GroupShuffleSplit`, `StandardScaler`
fit on train only, `penalty="l1"`, `solver="saga"`), per-module recovery-rate
analysis, and a final summary report table.

Key variables: `fold_change`, `mu_background`, `sigma2` (patient variance),
`dispersion` (NB overdispersion), `size_factor_mean/sd`, `C` (LASSO inverse
regularization strength), `n_modules`, `n_distractor_modules`, `genes_per_module`,
`within_module_correlation`. Random seed fixed at 42 throughout for reproducibility
(the copula draws use a separate `np.random.default_rng(seed)` so they don't
disturb the legacy global-state sequence used elsewhere).

## Correlated gene modules (step 2 — now implemented)
The 100 signal genes are split into 4 modules of 25 genes: 3 "disease modules"
(mean expression differs by group, same as the old flat design) + 1 "distractor
module" (just as internally correlated, but mean does NOT differ by group — a
negative control). Within-module correlation is induced via a **Gaussian copula**:
draw correlated standard-normal latents per (module, cell) using a compound-
symmetric correlation matrix, map to uniforms via the normal CDF, then invert each
gene's own NB CDF at that uniform (`nbinom.ppf`). This keeps every gene's marginal
NB distribution exactly as parameterized while inducing real pairwise correlation.
Background genes remain independent, unchanged.

**Validated result (current params, seed 42):** realized within-module correlation
(measured on size-factor-normalized expression, Group1 cells only, to strip the
technical size-factor floor — see gotchas below) ≈ 0.58-0.63 across disease modules
and the distractor module vs. target 0.6; background ≈ 0.20 (see "patient-effect
floor" gotcha below — this is expected, not leakage). Test accuracy 80.5% (train
93.4%) — in the target zone. **L1 recovers only 12-24% of each disease module's 25
genes** (M1: 24%, M2: 12%, M3: 20%) and **0% of the distractor module** — this is
the empirical evidence the thesis rests on: L1 doesn't recover whole correlated
modules, and it isn't fooled by correlation alone (only genuinely group-differential
modules get nonzero coefficients) — notably, this holds even though background and
distractor genes now carry real, non-trivial correlation (~0.2-0.63) from the
patient-effect floor, not just synthetic near-zero noise.

**Gotcha #1 — technical (size_factor) correlation floor:** raw expression counts
have a ~0.5 average pairwise correlation between essentially ALL genes (signal or
background), because `size_factor` multiplies every gene's mean within a cell — a
technical/sequencing-depth effect, not module structure. Always normalize by
`size_factor` (CPM-style) before measuring correlation to see the real
module-specific signal.

**Gotcha #2 — patient-effect correlation floor (added when `gamma_i` went onto
every gene):** `gamma_i` (the patient random effect) is applied additively to
*every* gene's mean in Step 8 — disease-module, distractor-module, and background
alike — since it represents a genuine patient-level effect that should shift all of
a patient's expression together, not just disease genes. Because `gamma_i` is
constant across all of a patient's cells, it survives the size_factor
normalization above and induces a real, nonzero correlation floor among ALL genes,
including background (~0.20 currently, exact value is a function of
`mu_background`/`sigma2`/`dispersion` together — read it fresh from
`background_corr` in the notebook, don't hand-copy a number into permanent docs).
**Read the distractor-vs-background contrast as "distractor sits meaningfully above
the shared floor (from the module copula); background sits at the floor" — not
"zero vs. correlated."** A `mu = max(mu, 1e-6)` floor guards against the additive
`gamma_i` form pushing a gene's mean non-positive for extreme draws (doesn't
trigger under current `sigma2`, but protects future larger-`sigma2` sweeps). Also
worth remembering for later: any future model that adds a patient covariate/random
effect will now soak up variance from every gene, not just disease genes.

## Parameter tuning notes (from experimentation)
- With the flat (pre-module) design: `fold_change=3, dispersion=15, sigma2=0.20,
  C=0.01` gave ~90% train/test accuracy.
- Increasing `fold_change`/`dispersion` and decreasing `sigma2` too far made the
  classification task nearly perfect (~99-100% accuracy) — too easy, not useful for
  showing L1's module-recovery failure.
- Target zone: ~75-85% test accuracy, with `C` high enough that coefficients aren't
  all zeroed out.
- With the module design (`fold_change=3, dispersion=15, sigma2=0.20,
  within_module_correlation=0.6`), `C=0.0054` lands at 80.5% test accuracy — retune
  `C` any time module/simulation params change (sweep e.g. 0.001 to 0.05; accuracy
  is very sensitive between C=0.005 and C=0.01 in this regime).

## Conventions
- Keep `np.random.seed(42)` for reproducibility unless intentionally testing across seeds.
- Long-format tables (`counts_df`) before pivoting wide (`expression_matrix`/`X_df`).
- Always fit `StandardScaler` on train split only, never on test or combined data.
