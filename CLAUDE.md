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

**Adjacent extension (not part of the 7 steps, but a natural follow-on):** disease
subtype/heterogeneity simulation — implemented in the notebook (see below) as
groundwork for a future contrastive-learning + clustering pass. Real diseases like
IPF/COPD/FHP usually aren't one uniform process; different patients sharing a
diagnosis can be driven by different (overlapping) gene programs. The notebook now
generates that ground truth (hidden per-patient subtype, subtype-gated module
boosting); actually building the contrastive-learning/clustering method to *recover*
subtypes from data is still future work, evaluated later with something like
Adjusted Rand Index / NMI against the hidden `subtype` column.

## Current state of `edits_syntheticdata.ipynb`
Simulator: imports/seed → parameters (incl. module + subtype knobs, all named
constants for previously-magic numbers) → patients_df (group + hidden `subtype`)
→ cells_df (size_factor) → genes_df (signal/background split + `module_id` +
`is_disease_module`) + `subtype_module_map` → gamma_i (patient random effect) →
group-specific means → **module-aware, subtype-gated Gaussian copula NB
simulation** (with correlation-matrix validity assertions and a mu-floor trigger
counter) → pivot to expression_matrix → summary stats → box plot → extended
data-QC/report section (dataset overview table, size-factor/gamma_i distributions,
expression histograms, per-module boxplots, per-module-by-subtype boxplots,
correlation heatmaps, PIT QQ plot, PCA). Followed by an L1 logistic regression
baseline test (train/test split by patient via `GroupShuffleSplit`, `StandardScaler`
fit on train only, `penalty="l1"`, `solver="saga"`, convergence check), per-module
recovery-rate analysis (with a disease-vs-distractor sanity assertion), a
true-effect-size-vs-coefficient recovery scatter plot, and a final summary report
table.

Key variables: `fold_change`, `mu_background`, `sigma2` (patient variance),
`dispersion` (NB overdispersion), `size_factor_mean/sd`, `n_modules`,
`n_distractor_modules`, `genes_per_module`, `within_module_correlation`,
`correlation_noise`, `n_subtypes`, `subtype_module_indices` (which disease-module
indices each subtype activates), `MU_FLOOR`/`EIGENVALUE_FLOOR`/`CHOLESKY_JITTER`/
`U_CLIP_EPS` (numerical-robustness constants), `L1_C`/`L1_TEST_SIZE`/`L1_MAX_ITER`.
All parameters live in one Step 2 config block. Random seed fixed at 42 throughout
for reproducibility (the copula draws and the subtype assignment each use a
separate `np.random.default_rng(seed)` so they don't disturb the legacy
global-state sequence used elsewhere). Notebook has no functions/module
extraction by design — stays a single linear notebook, everything runs top-to-bottom.

## Disease subtypes (heterogeneity ground truth, implemented 2026-07-23)
Each Group1 (diseased) patient secretly gets one of `n_subtypes` (default 3: S1,
S2, S3) subtype labels — ground truth only, never seen by any method below. Each
subtype activates a different, deliberately overlapping subset of the disease
modules (default: S1→{M1,M2}, S2→{M2,M3}, S3→{M1}). In Step 8, a disease module is
only boosted (fold_change applied) for a Group1 patient if their subtype activates
it; otherwise it sits at the same baseline as the distractor module. Group2 stays
homogeneous. Validated via a per-module-by-subtype boxplot (Step 11f-ii) — confirms
each subtype boosts exactly its assigned modules and nothing else.

**Consequence for accuracy/tuning:** since only a fraction of Group1 patients boost
any given module now, overall disease-vs-healthy signal is diluted and `C` needs to
move up (weaker regularization) relative to the pre-subtype design — retuned to
`C=0.0068` → 80.0% test accuracy (was `C=0.0054` pre-subtype). The C-sensitive
transition zone also shifted (now roughly 0.006-0.007, was 0.005-0.01) — expect to
re-sweep whenever subtype/module params change.

**Consequence for recovery rates:** a module tied to a rarer subtype gets less
training signal (fewer patients carry it) and is recovered less — e.g. M3 (S2 only,
3/10 Group1 patients) dropped to 0% recovery while M1/M2 (each ~7/10 patients) sit
around 32-36%. This is expected, not a bug: it directly shows that subtype rarity
makes standard L1's job harder, which is exactly the kind of realistic difficulty
the group-aware method (methodology step 4) should eventually help with.

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
(measured on size-factor-normalized expression, to strip the technical size-factor
floor — see gotchas below; disease modules measured within one activating
subtype's cells, distractor/background across all Group1 cells — see gotcha #3) ≈
0.57-0.63 across disease modules and the distractor module vs. target 0.6;
background ≈ 0.20 (see "patient-effect floor" gotcha below — this is expected, not
leakage). Test accuracy 80.0% (train 81.1%) — in the target zone, `C=0.0068`. **L1
recovers only 32-36% of M1/M2 and 0% of M3** (M3 is tied to the rarest subtype, S2,
only 3/10 Group1 patients — see "Disease subtypes" above) **and 0% of the
distractor module** — this is the empirical evidence the thesis rests on: L1
doesn't recover whole correlated modules, and it isn't fooled by correlation alone
(only genuinely group-differential modules get nonzero coefficients) — notably,
this holds even though background and distractor genes now carry real, non-trivial
correlation (~0.2-0.63) from the patient-effect floor, not just synthetic
near-zero noise.

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

**Gotcha #3 — subtype floor (added with disease-subtype heterogeneity):** for a
disease module, pooling ALL Group1 cells to measure correlation is now wrong —
only the subtypes that activate that module have it boosted, so pooling boosted
and unboosted patients together creates a large between-patient mean split (e.g.
M1 at ~boosted mean for S1/S3 patients, baseline for S2 patients) that inflates the
naive pooled correlation well past the copula's actual target (observed: ~0.87
instead of ~0.6 before this was fixed). Same failure mode as gotchas #1/#2, one
level up. Fix (already applied in the notebook): measure a disease module's
correlation only within cells from ONE subtype that activates it — every patient
in that subtype shares the same boost status for that module, so no between-patient
mean-split confound remains. Distractor/background are never subtype-boosted, so
pooling all of Group1 for them is still correct and unaffected.

## Parameter tuning notes (from experimentation)
- With the flat (pre-module, pre-subtype) design: `fold_change=3, dispersion=15,
  sigma2=0.20, C=0.01` gave ~90% train/test accuracy.
- Increasing `fold_change`/`dispersion` and decreasing `sigma2` too far made the
  classification task nearly perfect (~99-100% accuracy) — too easy, not useful for
  showing L1's module-recovery failure.
- Target zone: ~75-85% test accuracy, with `C` high enough that coefficients aren't
  all zeroed out.
- Module design, no subtypes (`fold_change=3, dispersion=15, sigma2=0.20,
  within_module_correlation=0.6`): `C=0.0054` → 80.5% test accuracy.
- Module + subtype design (same params, `n_subtypes=3`): signal is diluted since
  only a fraction of Group1 patients boost any given module, so `C` needed to move
  up — `C=0.0068` → 80.0% test accuracy. Retune `C` any time module/subtype/
  simulation params change; accuracy is very sensitive in a narrow band around
  whatever the current transition zone is (roughly 0.006-0.007 with subtypes,
  was 0.005-0.01 without) — always re-sweep rather than assuming the old zone holds.

## Conventions
- Keep `np.random.seed(42)` for reproducibility unless intentionally testing across seeds.
- Long-format tables (`counts_df`) before pivoting wide (`expression_matrix`/`X_df`).
- Always fit `StandardScaler` on train split only, never on test or combined data.
