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
**This has not been built yet — it is still the main remaining task** (methodology
step 4 below). Everything done so far is baseline-building and evidence-gathering
for why it's needed.

## Research questions
1. Do standard L1/stochastic gating fail to recover correlated gene modules?
2. Does incorporating correlation improve feature-selection stability/reproducibility?
3. Does group-aware gating improve interpretability without sacrificing accuracy?
4. Does it generalize across diseases/simulation settings?
5. (secondary) Does multi-view integration across cell types help patient representations?

## Methodology (7 steps)
1. Generate synthetic scRNA-seq data (Negative Binomial counts, patient/cell/gene
   hierarchy) — **done**.
2. Construct correlated gene modules — **done**.
3. Implement baselines (standard stochastic gates, L1, Elastic Net, + more) —
   **done, and expanded well beyond the original scope**: 10 feature-selection
   methods now implemented (see "Baseline Methods suite" below), evaluated
   per-module, not just signal-vs-background.
4. Build the group-aware gating framework — **NOT STARTED. This is the actual next
   task and the thesis's main contribution.** Nothing below this line exists yet.
5. Sweep simulation params (signal strength, correlation, noise, dimensionality) to
   stress-test robustness — not started (some ad hoc parameter tuning has happened,
   but not a systematic sweep).
6. Validate on real IPF/COPD/FHP datasets — not started.
7. Compare methods on AUROC/F1, feature recovery, stability, pathway enrichment —
   feature-recovery comparison done for the 10 baselines (see below); stability/
   pathway enrichment not started.

**Adjacent extension (not part of the 7 steps, but a natural follow-on):** disease
subtype/heterogeneity simulation — implemented as groundwork for a future
contrastive-learning + clustering pass. Real diseases like IPF/COPD/FHP usually
aren't one uniform process; different patients sharing a diagnosis can be driven by
different (overlapping) gene programs. The simulator generates that ground truth
(hidden per-patient subtype, subtype-gated module boosting); actually building the
contrastive-learning/clustering method to *recover* subtypes from data is still
future work, evaluated later with something like Adjusted Rand Index / NMI against
the hidden `subtype` column.

## Repo file map — READ THIS FIRST, the notebook situation is not simple
There are now **five notebooks** in this repo, representing two different lineages
that have diverged. If you're picking this project up fresh, do not assume any one
of these is "the" canonical file without checking — ask the user if unsure which
one they want worked on.

- **`edits_syntheticdata.ipynb`** — the original tracked notebook. Has the full
  module + subtype simulator (Gaussian copula, subtype-gated boosting) and an L1
  baseline with per-module/per-subtype recovery, a recovery scatter plot, and a
  summary report. **Has NOT received the TOC/header refactor or the 10-method
  Baseline Methods section** — it stopped evolving after commit `1b28e19` ("Add
  disease-subtype heterogeneity ground truth and notebook refactor"). Its own
  internal refactor (named constants, assertions, PIT QQ plot, etc.) predates and
  is separate from the later "refactor with comments" work described below.
- **`flat_syntheticdata.ipynb`** — no correlated modules, no subtypes. 100
  independent signal genes vs. 900 background. Has the TOC/header refactor.
  **As of 2026-07-30, this notebook is data-generation only** — the 10-method
  Baseline Methods section that used to live at the end of it was removed and
  moved into `baseline_methods_combined.ipynb` (see below and "Baseline Methods
  suite"). It now ends with a "Step 12: readiness check for `%run`" cell that
  asserts `counts_df`/`genes_df`/`patients_df`/`cells_df` exist and look
  well-formed, for whatever `%run`s it next. Generative params
  (`fold_change=3, dispersion=15, sigma2=0.20`) were manually updated by the user
  at some point to match the other two data notebooks, resolving an earlier
  parameter-mismatch concern (previously `fold_change=4, dispersion=50,
  sigma2=0.16`, which had made this dataset a much easier classification task than
  the other two, confounding cross-dataset comparisons).
- **`correlated_syntheticData.ipynb`** — note the **capital D** in "Data," this is
  the actual filename, not a typo. 4 gene modules (3 disease + 1 distractor), no
  subtypes. Has the TOC/header refactor. **As of 2026-07-30, data-generation only**
  (same restructuring as flat — Baseline Methods section removed, moved to
  `baseline_methods_combined.ipynb`, ends with the same readiness-check cell).
  Also has three untouched user-authored markdown cells at the very top
  (Date, Literature Review with 2 paper links, "Key bugs and their solutions") —
  **do not edit these**, the user explicitly asked they be preserved as-is.
  **Filesystem gotcha:** macOS/APFS is case-insensitive-but-case-preserving, so
  `correlated_syntheticdata.ipynb` (lowercase d) and `correlated_syntheticData.ipynb`
  (capital D) are THE SAME FILE on this machine — a save under one casing silently
  overwrites the other. This already destroyed an earlier STG section once. Always
  double check you're not creating an accidental case-only duplicate before writing
  a new file whose name differs from an existing one only by case.
- **`correlated_with_subtypes_syntheticdata.ipynb`** — 4 gene modules + per-patient
  disease subtypes. Has the TOC/header refactor. **As of 2026-07-30,
  data-generation only** (same restructuring — Baseline Methods section removed,
  moved to `baseline_methods_combined.ipynb`, ends with the same readiness-check
  cell). This is thematically the closest sibling to `edits_syntheticdata.ipynb`
  but is a **separately maintained, diverged copy** — changes to one do not
  propagate to the other.
- **`baseline_methods_combined.ipynb`** — **as of 2026-07-30, this is where all 10
  baseline methods' actual model code lives** (setup, training, evaluation) — it is
  no longer just a thin `%run`-and-checkpoint wrapper around embedded baseline
  sections in the data notebooks, since those sections were removed from the data
  notebooks and moved here. Structure, top to bottom:
  1. A pip-install cell (`!pip install -q numpy pandas scipy matplotlib seaborn
     scikit-learn xgboost torch stg group-lasso`) — added because Colab's default
     image lacks `stg` and `group-lasso`, which used to surface as a
     `ModuleNotFoundError` partway through a Colab run. No-op locally (everything
     already installed in `.venv`).
  2. "Colab Setup" section (clone `disease-heterogeneity` using a `GITHUB_TOKEN`
     Colab secret) — still a no-op locally, unchanged.
  3. One big `## Dataset N: ...` section per dataset (flat, correlated,
     correlated+subtypes), each containing, in order: a `%run` of that dataset's
     data notebook (with the local line active and the Colab-path line commented
     out, swap which is commented when porting to Colab) → a "Step 12: setup" cell
     (builds `X_df`/train-test split, same as the old embedded setup cells used to)
     → ten `### Method N: ...` mini-sections (one per baseline method) → a
     "Step 23: full comparison" cell → a checkpoint cell that saves that dataset's
     `baseline_df`/`module_recovery_df` to CSV → `%reset -f` to fully clear state
     before the next dataset (prevents e.g. `module_names` leaking across
     datasets).
  4. A final "Combined comparison across all three datasets" section that reloads
     the checkpoint CSVs and produces the cross-dataset comparison (table + bar
     charts).

## Collaborative workflow (established 2026-07-29)
The user develops/edits in **Google Colab**, connected to a collaborator's GitHub
repo: `carmen666666999/disease-heterogeneity`, branch `Sierra` (this is shared
with "the Yale people" — treat it as at least semi-private/collaborative, not the
user's own solo repo). The loop:

1. User hands Claude a current copy of the notebook(s) (or Claude works from what's
   already in this local repo).
2. Claude edits locally in `sierra-gating-project/`, verifies with a real execution
   (`jupyter nbconvert --to notebook --execute`), commits and pushes to
   `sierramendoza543/sierra-gating-project` on GitHub.
3. User reviews the diff via the **ReviewNB** GitHub App (gives a rendered,
   cell-by-cell diff instead of raw JSON — this is *why* commits go to
   `sierra-gating-project` at all instead of directly to the Colab-linked repo).
4. User manually re-applies the changes they like into the real Colab notebook
   (which lives in/against `carmen666666999/disease-heterogeneity`).

**Practical consequences of this loop:**
- Local commits to `sierra-gating-project` are a staging/review area, not
  necessarily the final home of this code — don't assume "committed here" means
  "done," it means "ready for the user to review and port."
- `%run` (Colab/IPython magic) needs an actual local file path — opening a notebook
  from GitHub in Colab does NOT place a corresponding file on the Colab VM's disk.
  The established pattern for the real Colab notebooks: clone (or download) from
  GitHub at the top of the notebook, then `%run` the local cloned path.
- `disease-heterogeneity` is not the user's own repo, so a **fine-grained** GitHub
  PAT can't be scoped to it (the "Resource owner" dropdown only offers your own
  account + orgs you belong to, never another individual's personal account). The
  working solution is a **classic PAT** with the `repo` scope, stored in Colab's
  Secrets manager (key icon 🔑, name `GITHUB_TOKEN`), retrieved via
  `google.colab.userdata.get('GITHUB_TOKEN')`. Needs per-notebook "Notebook access"
  granted once for each new Colab notebook that wants to read it.
- Every dataset the Colab-side automation needs, `%run` needs its own line pointing
  at the cloned path (e.g. `%run disease-heterogeneity/flat_syntheticdata.ipynb`) —
  this is already stubbed out (commented) in `baseline_methods_combined.ipynb`.

## Baseline Methods suite (10 methods, added 2026-07-29)
All 10 methods' setup/training/evaluation code lives in
`baseline_methods_combined.ipynb` (see "Repo file map" above for the exact
structure) — one big section per dataset, each starting with a "Step 12: setup"
cell that rebuilds `X_df`/train-test split from that dataset's `counts_df`/
`genes_df`/`cells_df`, followed by a mini-section per method. **This moved out of
the three data notebooks on 2026-07-30** — if you're looking at an older mental
model (or an older commit) where the data notebooks each had their own embedded
Baseline Methods section, that's stale; don't assume it without checking. Methods,
grouped by paradigm:

- **Embedded-linear**: L1 logistic regression (the original headline baseline),
  ElasticNet, Ridge, Group Lasso, L1-penalized linear SVM.
- **Embedded-nonlinear**: Random Forest, XGBoost.
- **Wrapper**: Recursive Feature Elimination (RFE).
- **Filter**: mutual-information univariate selection (`SelectKBest`).
- **Gate-based**: STG (Stochastic Gates, Yamada et al. ICML 2020, via the
  official `stg` pip package).

**How "selected"/"recovered" is defined, uniformly across all 10:** methods with
genuine built-in sparsity (L1, ElasticNet, Group Lasso, L1-SVM, STG) use their own
nonzero-coefficient / gate-open-probability-over-0.5 threshold. Methods with no
natural sparsity (Ridge, Random Forest, XGBoost) or that need a fixed budget upfront
(RFE, SelectKBest) instead take the **top-100** genes by importance score (100 =
the true number of signal-carrying genes in the simulation). **Caveat not yet
written into the notebooks:** this means the top-100 methods are being handed a
much more generous selection budget than L1/STG's self-chosen sparsity (typically
10-40 genes) — so "top-100 methods recover modules better" partly reflects budget
size, not purely method quality. Worth a one-line clarifying note in each
notebook's final comparison cell next time this section is touched.

**Group Lasso is the most thesis-relevant baseline**: it's an existing, simple form
of "group-aware" sparsity (Yuan & Lin 2006) — exactly the kind of thing the
eventual group-aware gating framework (methodology step 4) needs to be compared
against. On the two module-based notebooks it uses the real `module_id` groups
(background genes each get their own singleton group); on the flat notebook it
necessarily reduces to per-gene behavior since there are no real groups there.

**Headline result (correlated notebook, 10-method comparison):** L1 and STG show
the LOWEST per-module recovery of all 10 methods (roughly 10-24%), vs. 75-100% for
the top-100-budget methods (RF/RFE/Ridge/SelectKBest) — and every single method
gives ~0% recovery to the distractor module. This extends the thesis's core
evidence from "L1 alone fails to recover modules" to "a whole battery of standard
methods fails," which is a stronger empirical foundation for methodology step 4.
On the subtypes notebook, STG's behavior changed qualitatively under signal
dilution (much less sparse, higher recovery, at the identical `lam` that gave very
sparse/low-recovery behavior on the plain correlated notebook) — a real, reported
finding, not something to average away.

**IMPORTANT interpretation correction — do not conflate these two things when
writing up or re-deriving "why baselines fail":**
- **The distractor module's ~0% recovery across all 10 methods is NOT evidence of
  failure.** M4 has zero true classification signal by design (it's a negative
  control — just as correlated as a real module, but mean does not differ by
  group). A supervised method correctly assigning it ~0% recovery is *correct,
  expected behavior*, not a shortcoming to fix. This was stated as "failure
  evidence" earlier in this project's analysis and that framing was wrong — it was
  corrected once already; don't reintroduce it.
- **The actual failure is under-recovery of the REAL disease modules (M1-M3) by
  L1/STG/L1-SVM at their own self-chosen sparsity.** These methods correctly
  determine a module matters, then arbitrarily keep only 3-6 of a module's 25
  genes and drop the rest — that's the literal "picks 1-2 representative genes,
  drops the rest" problem from this project's core problem statement, and it's a
  real, measured result (not an artifact). This is the evidence methodology step 4
  needs to beat, not the distractor-module numbers.

**Budget-artifact follow-up test (SelectKBest-MI, tested 2026-07-30, not yet
folded into any notebook):** the top-100 budget handed to Ridge/RF/XGBoost/RFE/
SelectKBest-MI is itself a form of oracle knowledge (100 = the true number of
signal-carrying genes, which you'd never know in a real dataset) — Ridge and
Random Forest in particular have *no* sparsity mechanism at all without an
external cutoff (L2 coefficients and Gini/gain importances are never exactly
zero), so their "100% recovery" is entirely an artifact of the budget chosen for
them, not evidence of understanding modules. To check whether SelectKBest-MI
specifically has a legitimate non-oracle alternative, it was re-run with no fixed
`k` at all: instead, a permutation-based significance threshold (shuffle the
labels 100 times, take the 95th/99th percentile of the resulting null max-MI
score across all 1000 genes — a standard max-statistic FWER-controlled test).
Result, on **both** the correlated and correlated+subtypes datasets: it
self-selected **75 genes** with **100% recovery on M1, M2, and M3, 0% on the
distractor, and 0% background leak** — cleaner than the budgeted version (which
leaked 8% distractor + 2.56% background because it was forced to fill 100 slots),
and, notably, still 100% on M3 in the subtypes dataset (the module diluted to only
3/10 Group1 patients, where L1 collapses to 0%). Mechanistic reason it survives
dilution where L1 doesn't: MI scores each gene in total isolation, so it never
competes against other genes for a shared sparsity budget the way L1's
coefficients do — a rarer signal doesn't get crowded out. **Practical implication:
SelectKBest-MI should NOT be grouped with Ridge/RF as "only works because of an
oracle budget"** — unlike them, it has a real, standard, budget-free alternative
and it holds up. RFE's non-oracle alternative would be RFECV (stops once
cross-validated accuracy plateaus); this was *not* tested, but is predicted to
collapse early and look more like L1's failure, since L1 already hits 93-99%
accuracy using only 10-35 genes, well before all 75 true genes are needed for
accuracy alone. Code and full output for the SelectKBest-MI test: standalone
script (not a notebook cell) plus a write-up in `SelectKBest-MI_permutation_test.docx`
in the repo root — **this file is not committed to git or tracked, ask the user
if you need its content and it's not present.**

**Per-notebook tuned hyperparameters** (retune if Step 2's generative params
change): L1 `C` — flat 0.01, correlated 0.0054, subtypes 0.0068. STG `lam` — flat
0.5, correlated/subtypes 0.2. ElasticNet/LinearSVM-L1 reuse the same `C` as that
notebook's L1 for simplicity (not independently tuned). Group Lasso: `group_reg=0.05,
l1_reg=0.0` on all three. Random Forest/XGBoost: `n_estimators=300`. RFE:
`step=50` (NOT `step=1` — with 1000 features, step=1 means ~1000 refits and is
extremely slow). All three use the same `GroupShuffleSplit(test_size=0.2,
random_state=42)` split-by-patient as the data-generation-adjacent L1 test that
predates this section.

## Known bugs / environment gotchas
- **PyTorch + XGBoost segfault (macOS, likely other platforms too).** Importing
  `torch` before `xgboost` in the same process, then actually running both, causes
  a silent SIGSEGV (exit code 139, no traceback — looks like the process just
  vanished). Both bundle their own OpenMP runtime and conflict once both execute.
  Fix, baked into every "Baseline Methods setup" cell: set
  `os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"` and `os.environ["OMP_NUM_THREADS"]="1"`
  **before any other imports**, import `xgboost` before `torch`, call
  `torch.set_num_threads(1)`, and pass `n_jobs=1` to XGBoost/RandomForest. Locally
  this machine also needed `brew install libomp` (XGBoost's dylib dependency) —
  Colab ships with this preinstalled, so that specific step is Mac-local-only, not
  a notebook-code issue.
- **`stg` package + Python 3.10+.** The package (last released a while ago) still
  references `collections.Sequence`/`Mapping`/`Iterable`/`Set`, removed from
  `collections` in Python 3.10 (moved to `collections.abc`). Needs a compatibility
  shim at import time in every environment, including Colab:
  ```python
  import collections, collections.abc
  for _name in ("Iterable", "Mapping", "Sequence", "Set", "MutableMapping",
                "MutableSequence", "MutableSet", "Callable", "Hashable",
                "Sized", "Container"):
      if not hasattr(collections, _name):
          setattr(collections, _name, getattr(collections.abc, _name))
  ```
- **`STG`'s `random_state` constructor arg is accepted but never actually used
  internally** (checked the installed package source) — it relies on the global
  PyTorch/NumPy RNG state instead. Call `torch.manual_seed(42)` and
  `np.random.seed(42)` explicitly right before `.fit()` for reproducibility.
- **`STG.evaluate()` only prints a summary and returns `None`** — compute accuracy
  from `.predict()` manually: `(model.predict(X) == y).mean()`.
- **scikit-learn ≥1.8 deprecation noise**: emits a `FutureWarning` about `penalty=`
  being replaced by `l1_ratio=`. Harmless for this project's pinned
  `penalty="l1"/solver="saga"` usage; suppressed via
  `warnings.filterwarnings("ignore", category=FutureWarning)` in the relevant setup
  cells so it doesn't clutter output (and doesn't leak local file paths into
  committed notebook outputs via the warning's traceback).

## Correlated gene modules — design and correlation-measurement gotchas
The 100 signal genes are split into 4 modules of 25 genes: 3 "disease modules"
(mean expression differs by group) + 1 "distractor module" (just as internally
correlated, but mean does NOT differ by group — a negative control). Within-module
correlation is induced via a **Gaussian copula**: draw correlated standard-normal
latents per (module, cell) using a compound-symmetric correlation matrix
(`within_module_correlation`, optionally jittered by `correlation_noise` — see
`build_module_corr_matrix`), map to uniforms via the normal CDF, then invert each
gene's own NB CDF at that uniform (`nbinom.ppf`). This keeps every gene's marginal
NB distribution exactly as parameterized while inducing real pairwise correlation.
Background genes remain independent.

Three correlation-measurement gotchas were found, each the same failure mode one
level up from the last — **the general lesson: any new axis of patient
heterogeneity tends to create its own correlation-measurement confound the first
time it's added; check for this before trusting a new "realized correlation" number**:

1. **Technical (size_factor) floor**: raw counts have a ~0.5 average pairwise
   correlation between essentially ALL genes, because `size_factor` multiplies
   every gene's mean within a cell. Fix: normalize by `size_factor` (CPM-style)
   before measuring correlation.
2. **Patient-effect floor**: `gamma_i` (patient random effect) is additive on
   *every* gene (disease, distractor, background alike, since it represents a
   genuine patient-level effect). It's constant within a patient, survives the
   size_factor normalization, and induces its own nonzero correlation floor across
   ALL genes (~0.20 on background, exact value depends on
   `mu_background`/`sigma2`/`dispersion` together — read it fresh from
   `background_corr`, don't hand-copy a number anywhere permanent). Read
   distractor-vs-background as "distractor sits meaningfully above the shared
   floor; background sits at the floor," not "zero vs. correlated." A
   `mu = max(mu, MU_FLOOR)` guards against this additive form pushing a mean
   non-positive (doesn't trigger under current `sigma2`, protects future sweeps).
3. **Subtype floor** (subtypes notebooks only): pooling ALL Group1 cells to
   measure a disease module's correlation is wrong once subtypes exist — only the
   subtypes that activate that module have it boosted, so pooling boosted/unboosted
   patients together inflates the naive correlation (~0.87 observed vs. 0.6
   target). Fix: measure a disease module's correlation only within cells from ONE
   subtype that activates it (every patient in that subtype shares the same boost
   status, so no between-patient mean-split confound). Distractor/background are
   never subtype-boosted, so pooling all of Group1 is still correct for them.

## Disease subtypes (heterogeneity ground truth)
Each Group1 (diseased) patient secretly gets one of `n_subtypes` (default 3: S1,
S2, S3) subtype labels — ground truth only, never seen by any method. Each subtype
activates a different, deliberately overlapping subset of the disease modules
(default: S1→{M1,M2}, S2→{M2,M3}, S3→{M1}). A disease module is only boosted
(fold_change applied) for a Group1 patient if their subtype activates it;
otherwise it sits at the same baseline as the distractor module. Group2 stays
homogeneous. A module tied to a rarer subtype gets less training signal and is
recovered less by standard methods (e.g. M3/S2, only 3/10 Group1 patients, sits
near 0% L1 recovery) — expected/realistic, not a bug, and exactly the kind of
difficulty the eventual group-aware method should help with.

## Parameter tuning notes
- **Generative params are now matched across all three data notebooks**:
  `fold_change=3, dispersion=15, sigma2=0.20, mu_background=2,
  size_factor_mean=2000, size_factor_sd=300, n_patients=10/10,
  cells_per_patient=50, n_genes_total=1000, n_signal_genes=100`. This was a
  deliberate fix (flat used to run easier params) so cross-dataset comparisons
  aren't confounded by task difficulty — only the intended structural variable
  (modules present/absent, subtypes present/absent) should differ between them.
- **Module-specific params** (correlated + subtypes notebooks): `n_modules=4,
  n_distractor_modules=1, within_module_correlation=0.6, correlation_noise=0.0`.
- **Target zone for L1's own accuracy**: ~75-85% test accuracy is the useful
  difficulty level for showing module-recovery failure; ~99-100% is "too easy" and
  shows nothing interesting. See "Baseline Methods suite" above for the other 9
  methods' hyperparameters, which are NOT all individually tuned to this same zone
  (only L1 and STG were carefully swept; the rest use standard defaults).
- Retune `C`/`lam`/etc. any time Step 2's generative parameters change —
  accuracy is sensitive in a narrow band around whatever the current transition
  zone is; always re-sweep rather than assuming an old tuning holds.

## Conventions
- Keep `np.random.seed(42)` for reproducibility unless intentionally testing across
  seeds. Copula draws and subtype assignment each use a separate
  `np.random.default_rng(seed)` so they don't disturb the legacy global-state
  sequence used elsewhere (size_factor, gamma_i, etc.).
- Long-format tables (`counts_df`) before pivoting wide (`expression_matrix`/`X_df`).
- Always fit `StandardScaler` on train split only, never on test or combined data.
- **No functions, no `.py` module extraction** — explicit prior user direction.
  These stay single, linear notebooks that run top-to-bottom; don't refactor logic
  into importable helpers even if it would reduce duplication across the three
  data notebooks.
- **Documentation style** (established in the "refactor with comments" pass): a
  Table of Contents cell near the top; every step gets a real markdown header
  (`## Step N: ...`, not a bare "Step N" one-liner) so Colab's/GitHub's
  auto-generated outline picks it up; a "quick tuning guide" + debugging notes near
  Step 2 (parameters) and Step 8 (the NB simulation); comments explain *why*, not
  just *what*. `flat_syntheticdata.ipynb`/`correlated_syntheticData.ipynb`/
  `correlated_with_subtypes_syntheticdata.ipynb` all follow this now;
  `edits_syntheticdata.ipynb` does not (see "Repo file map" above).
- Before creating any new file, check whether its name collides case-insensitively
  with an existing one (see the filesystem gotcha above) — `ls` the directory
  first if in doubt.
- Each of the three data notebooks ends with a "Step 12: readiness check for
  `%run`" cell (added 2026-07-30) — asserts `counts_df`/`genes_df`/`patients_df`/
  `cells_df` exist and look well-formed (right columns, no NaNs, unique ids), with
  conditional checks for `module_names`/`disease_modules` and `subtype` that only
  fire on the notebooks where those exist. Identical cell content works unmodified
  across all three notebooks — don't hand-write a different version per notebook.
- Local Python environment: a `.venv/` in the repo root (gitignored), with
  `numpy pandas scipy matplotlib seaborn scikit-learn jupyter nbconvert stg torch
  xgboost group-lasso` installed. Verify any notebook change with a real
  `jupyter nbconvert --to notebook --execute` run before committing — this project
  has a strong pattern of "looks right on paper, fails at runtime" issues (the
  OpenMP segfault, the `collections` shim, `%run`'s cwd-follows-notebook-location
  behavior) that only show up on actual execution.

## What's actually next
1. **Methodology step 4 — the group-aware gating framework itself.** This is the
   whole point of the thesis and hasn't been started. Likely direction based on
   the framing above: a stochastic-gates-style architecture where gates are
   coupled/shared within a module (known or inferred) rather than independent
   per-gene, trained jointly with the classifier. The evidence this needs to beat
   is L1/STG/L1-SVM's low *self-chosen-sparsity* recovery of the real disease
   modules (M1-M3) — NOT the distractor module's 0% (see the "IMPORTANT
   interpretation correction" note above, which flags a mistake already made once
   in this project's own analysis — don't repeat it).
2. Add the budget-asymmetry caveat (top-100 vs. self-chosen sparsity) to the
   Baseline Methods comparison cells — still not done in the notebook itself. The
   SelectKBest-MI permutation-threshold variant (see above) has now been tested
   standalone and could be folded into `baseline_methods_combined.ipynb` as an
   eleventh method/variant if the user wants it there rather than just in the
   separate `.docx` write-up.
3. Decide what to do about `edits_syntheticdata.ipynb` vs.
   `correlated_with_subtypes_syntheticdata.ipynb` having diverged — ask the user
   rather than assuming one supersedes the other.
4. Methodology step 5 (systematic parameter sweeps) and step 6 (real IPF/COPD/FHP
   data validation) haven't been started at all. One specific, well-motivated sweep
   direction raised by the user (2026-07-30): raising `within_module_correlation`
   (currently 0.6) rather than generically increasing noise/difficulty — the
   reasoning being that higher intra-module correlation should make L1-style
   "arbitrarily keep one representative gene" pruning worse, while ElasticNet/Group
   Lasso (built to keep correlated genes together) should stay robust or improve,
   sharpening the specific contrast the thesis is trying to show rather than just
   making the task harder for everyone uniformly. Not yet run.
5. **Secondary/adjacent scoping question, resolved 2026-07-30**: the bigger Yale
   lab project integrates across cell types (research question 5); this smaller
   thesis project's simulator has no cell-type dimension at all and, per user
   decision, should NOT take on cell-type integration now — it would be a
   substantial separate methodological addition (typed cells, per-cell-type
   modules, a multi-view architecture) that would dilute focus away from step 4,
   which hasn't been built yet. Treat as future work, not current scope, unless
   the user revisits this.
