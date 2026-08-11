"""
data_generators.py

Standalone, parameterized data-generation functions for the synthetic
scRNA-seq framework (Group-Aware Stochastic Gating project).

Two functions:
    generate_flat_data(...)         -> uncorrelated (flat) synthetic dataset
    generate_correlated_data(...)   -> module-correlated synthetic dataset
                                        (Gaussian-copula NB sampling)

Both return a dict containing the same DataFrames the original notebooks
produce (counts_df, genes_df, patients_df, cells_df), plus a ready-to-use
wide expression matrix X (cells x genes) and label vector y
(0 = Group2/control, 1 = Group1/patient), for direct use with any
scikit-learn-style model, L1/LASSO, or stochastic gating.

Usage:
    from data_generators import generate_flat_data, generate_correlated_data

    flat = generate_flat_data()

    corr_02 = generate_correlated_data(
        n_modules=1, n_distractor_modules=0, within_module_correlation=0.2
    )
"""

import numpy as np
import pandas as pd
from scipy.stats import nbinom, norm


def generate_flat_data(
    n_patients_group1=10,
    n_patients_group2=10,
    cells_per_patient=50,
    n_genes_total=1000,
    n_signal_genes=100,
    fold_change=3,
    mu_background=2,
    sigma2=0.20,
    dispersion=15.0,
    size_factor_mean=2000,
    size_factor_sd=300,
    seed=42,
):
    """
    Regenerate the flat (uncorrelated) synthetic scRNA-seq dataset.

    Returns
    -------
    dict with keys:
        counts_df, genes_df, patients_df, cells_df : pandas DataFrames,
            same shape/columns as in flat_syntheticdata.ipynb
        X : wide cells x genes expression matrix (pandas DataFrame)
        y : per-cell group label aligned to X.index
            (0 = Group2/control, 1 = Group1/patient)
    """
    rng = np.random.default_rng(seed)

    # Patients
    patient_ids = [f"P{i+1}" for i in range(n_patients_group1 + n_patients_group2)]
    groups = ["Group1"] * n_patients_group1 + ["Group2"] * n_patients_group2
    patients_df = pd.DataFrame({"patient_id": patient_ids, "group": groups})
    patients_df["gamma_i"] = rng.normal(0, np.sqrt(sigma2), len(patients_df))

    # Cells / size factors
    cell_rows = []
    for _, row in patients_df.iterrows():
        size_factors = np.clip(
            rng.normal(size_factor_mean, size_factor_sd, cells_per_patient), 1, None
        )
        for c_idx in range(cells_per_patient):
            cell_rows.append({
                "patient_id": row.patient_id,
                "group": row.group,
                "cell_id": f"{row.patient_id}_C{c_idx+1}",
                "size_factor": size_factors[c_idx],
            })
    cells_df = pd.DataFrame(cell_rows)

    # Genes
    gene_ids = [f"G{i+1}" for i in range(n_genes_total)]
    signal_gene_ids = gene_ids[:n_signal_genes]
    background_gene_ids = gene_ids[n_signal_genes:]
    genes_df = pd.DataFrame({
        "gene_id": gene_ids,
        "gene_type": ["signal"] * n_signal_genes
        + ["background"] * (n_genes_total - n_signal_genes),
    })

    # Group-specific means
    mu_group1_signal = mu_background * fold_change
    mu_group2_signal = mu_background

    # Simulate counts (gamma_i applied to every gene, signal and background)
    count_rows = []
    for _, cell in cells_df.iterrows():
        gamma_i = patients_df.loc[
            patients_df.patient_id == cell.patient_id, "gamma_i"
        ].values[0]
        signal_mu = (
            mu_group1_signal if cell.group == "Group1" else mu_group2_signal
        ) + gamma_i
        background_mu = mu_background + gamma_i

        for gene_id in signal_gene_ids:
            mu = cell.size_factor * signal_mu
            var = mu + mu**2 / dispersion
            p = mu / var
            n = mu * p / (1 - p)
            count_rows.append({
                **cell.to_dict(),
                "gene_id": gene_id,
                "gene_type": "signal",
                "expression": nbinom.rvs(n, p, random_state=rng),
            })
        for gene_id in background_gene_ids:
            mu = cell.size_factor * background_mu
            var = mu + mu**2 / dispersion
            p = mu / var
            n = mu * p / (1 - p)
            count_rows.append({
                **cell.to_dict(),
                "gene_id": gene_id,
                "gene_type": "background",
                "expression": nbinom.rvs(n, p, random_state=rng),
            })

    counts_df = pd.DataFrame(count_rows)

    X = counts_df.pivot_table(index="cell_id", columns="gene_id", values="expression")
    y = cells_df.set_index("cell_id").loc[X.index, "group"].map(
        {"Group2": 0, "Group1": 1}
    )

    return {
        "counts_df": counts_df,
        "genes_df": genes_df,
        "patients_df": patients_df,
        "cells_df": cells_df,
        "X": X,
        "y": y,
    }


def generate_correlated_data(
    n_patients_group1=10,
    n_patients_group2=10,
    cells_per_patient=50,
    n_genes_total=1000,
    n_signal_genes=100,
    fold_change=3,
    mu_background=2,
    sigma2=0.20,
    dispersion=15.0,
    size_factor_mean=2000,
    size_factor_sd=300,
    n_modules=1,                    # simple scenario default: 1 module
    n_distractor_modules=0,         # simple scenario default: no distractor
    within_module_correlation=0.6,  # <-- the parameter to sweep (0.2/0.4/0.6/0.8)
    correlation_noise=0.0,
    seed=42,
):
    """
    Regenerate the module-correlated synthetic scRNA-seq dataset
    (Gaussian-copula NB sampling).

    Returns
    -------
    dict with the same keys as generate_flat_data(), plus:
        module_names, disease_modules, distractor_module_names
    """
    assert n_signal_genes % n_modules == 0, (
        "n_signal_genes must divide evenly by n_modules"
    )
    genes_per_module = n_signal_genes // n_modules
    MU_FLOOR = 1e-6
    rng = np.random.default_rng(seed)

    # Patients / cells -- identical mechanism to generate_flat_data
    patient_ids = [f"P{i+1}" for i in range(n_patients_group1 + n_patients_group2)]
    groups = ["Group1"] * n_patients_group1 + ["Group2"] * n_patients_group2
    patients_df = pd.DataFrame({"patient_id": patient_ids, "group": groups})
    patients_df["gamma_i"] = rng.normal(0, np.sqrt(sigma2), len(patients_df))

    cell_rows = []
    for _, row in patients_df.iterrows():
        size_factors = np.clip(
            rng.normal(size_factor_mean, size_factor_sd, cells_per_patient), 1, None
        )
        for c_idx in range(cells_per_patient):
            cell_rows.append({
                "patient_id": row.patient_id,
                "group": row.group,
                "cell_id": f"{row.patient_id}_C{c_idx+1}",
                "size_factor": size_factors[c_idx],
            })
    cells_df = pd.DataFrame(cell_rows)

    # Genes / modules
    gene_ids = [f"G{i+1}" for i in range(n_genes_total)]
    signal_gene_ids = gene_ids[:n_signal_genes]
    background_gene_ids = gene_ids[n_signal_genes:]
    module_names = [f"M{i+1}" for i in range(n_modules)]
    disease_modules = module_names[: n_modules - n_distractor_modules]
    distractor_module_names = module_names[n_modules - n_distractor_modules:]

    gene_module = {}
    for m_idx, module in enumerate(module_names):
        start = m_idx * genes_per_module
        for gene_id in signal_gene_ids[start:start + genes_per_module]:
            gene_module[gene_id] = module

    genes_df = pd.DataFrame({
        "gene_id": gene_ids,
        "gene_type": ["signal"] * n_signal_genes
        + ["background"] * (n_genes_total - n_signal_genes),
    })
    genes_df["module_id"] = genes_df.gene_id.map(gene_module).fillna("none")
    module_gene_lists = {
        m: genes_df.loc[genes_df.module_id == m, "gene_id"].tolist()
        for m in module_names
    }

    # Module correlation matrix + Cholesky factor (shared across modules)
    corr = np.eye(genes_per_module) * (1 - within_module_correlation) + within_module_correlation
    if correlation_noise > 0:
        jitter = rng.normal(0, correlation_noise, (genes_per_module, genes_per_module))
        jitter = (jitter + jitter.T) / 2
        np.fill_diagonal(jitter, 0)
        corr = corr + jitter
        eigvals, eigvecs = np.linalg.eigh(corr)
        eigvals = np.clip(eigvals, 1e-6, None)
        corr = eigvecs @ np.diag(eigvals) @ eigvecs.T
        d = np.sqrt(np.diag(corr))
        corr = corr / np.outer(d, d)
    corr = corr + np.eye(genes_per_module) * 1e-8
    chol = np.linalg.cholesky(corr)

    mu_group1_signal = mu_background * fold_change
    mu_group2_signal = mu_background

    count_rows = []
    for _, cell in cells_df.iterrows():
        gamma_i = patients_df.loc[
            patients_df.patient_id == cell.patient_id, "gamma_i"
        ].values[0]
        disease_mu_base = mu_group1_signal if cell.group == "Group1" else mu_group2_signal

        for module in module_names:
            base = disease_mu_base if module in disease_modules else mu_background
            mu_gene = max(cell.size_factor * (base + gamma_i), MU_FLOOR)
            var = mu_gene + mu_gene**2 / dispersion
            p = mu_gene / var
            n = mu_gene * p / (1 - p)

            z = chol @ rng.standard_normal(genes_per_module)
            u = np.clip(norm.cdf(z), 1e-6, 1 - 1e-6)
            exprs = nbinom.ppf(u, n, p).astype(int)

            for gene_id, expr in zip(module_gene_lists[module], exprs):
                count_rows.append({
                    **cell.to_dict(),
                    "gene_id": gene_id,
                    "gene_type": "signal",
                    "module_id": module,
                    "expression": expr,
                })

        for gene_id in background_gene_ids:
            mu = max(cell.size_factor * (mu_background + gamma_i), MU_FLOOR)
            var = mu + mu**2 / dispersion
            p = mu / var
            n = mu * p / (1 - p)
            count_rows.append({
                **cell.to_dict(),
                "gene_id": gene_id,
                "gene_type": "background",
                "module_id": "none",
                "expression": nbinom.rvs(n, p, random_state=rng),
            })

    counts_df = pd.DataFrame(count_rows)

    X = counts_df.pivot_table(index="cell_id", columns="gene_id", values="expression")
    y = cells_df.set_index("cell_id").loc[X.index, "group"].map(
        {"Group2": 0, "Group1": 1}
    )

    return {
        "counts_df": counts_df,
        "genes_df": genes_df,
        "patients_df": patients_df,
        "cells_df": cells_df,
        "module_names": module_names,
        "disease_modules": disease_modules,
        "distractor_module_names": distractor_module_names,
        "X": X,
        "y": y,
    }
