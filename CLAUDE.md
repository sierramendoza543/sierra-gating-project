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
2. Construct correlated gene modules — **NOT done yet, see "Current gap" below.**
3. Implement baselines: standard stochastic gates, L1, Elastic Net — **L1 done**.
4. Build the group-aware gating framework — **not started**, the main contribution.
5. Sweep simulation params (signal strength, correlation, noise, dimensionality) to
   stress-test robustness.
6. Validate on real IPF/COPD/FHP datasets.
7. Compare methods on AUROC/F1, feature recovery, stability, pathway enrichment.

## Current state of `edits_syntheticdata.ipynb`
11-step simulator: imports/seed → parameters → patients_df → cells_df (size_factor)
→ genes_df (signal/background split) → gamma_i (patient random effect) →
group-specific means (fold_change applied here) → NB count simulation → pivot to
expression_matrix → summary stats → box plot. Followed by an L1 logistic regression
baseline test (train/test split, StandardScaler fit on train only, `penalty="l1"`,
`solver="saga"`).

Key variables: `fold_change`, `mu_background`, `sigma2` (patient variance),
`dispersion` (NB overdispersion), `size_factor_mean/sd`, `C` (LASSO inverse
regularization strength). Random seed fixed at 42 throughout for reproducibility.

## Current gap (the next real task)
All 100 "signal" genes currently share one flat structure — there are no distinct
correlated **modules** within them yet. Need to split signal genes into several
groups (e.g. 4 modules × 25 genes) where genes *within* a module share a common
latent factor (correlated), but modules are independent of each other. Then rerun
the L1 baseline **per module** to check whether it recovers a whole module or just
1-2 representative genes — this is the empirical evidence the whole thesis rests on.

## Parameter tuning notes (from experimentation)
- Original settings (`fold_change=3, dispersion=15, sigma2=0.20, C=0.01`) give
  ~90% train/test accuracy — a reasonable baseline difficulty.
- Increasing `fold_change`/`dispersion` and decreasing `sigma2` too far made the
  classification task nearly perfect (~99-100% accuracy) — too easy, not useful for
  showing L1's module-recovery failure.
- Target zone: ~75-85% test accuracy, with `C` high enough that coefficients aren't
  all zeroed out. Sweep `C` (e.g. 0.005 to 1.0) whenever data parameters change.

## Conventions
- Keep `np.random.seed(42)` for reproducibility unless intentionally testing across seeds.
- Long-format tables (`counts_df`) before pivoting wide (`expression_matrix`/`X_df`).
- Always fit `StandardScaler` on train split only, never on test or combined data.
