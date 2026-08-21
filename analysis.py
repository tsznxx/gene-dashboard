# analysis.py

# for test
import streamlit as st

import numpy as np
import pandas as pd
from pycombat import Combat
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==================================================
# Batch Correction
# ==================================================


def apply_combat(expr_df, meta_df, batch_column):
    """
    Apply ComBat batch correction.

    Parameters
    ----------
    expr_df : pd.DataFrame

        Format:
        Gene    Sample1 Sample2 ...

    meta_df : pd.DataFrame

        Must contain:
        Sample
        batch_column

    batch_column : str

        Metadata column defining batch.

    Returns
    -------
    pd.DataFrame
        Batch-corrected expression matrix
        in the same format as input.
    """

    #
    # Expression matrix
    #
    expr_mat = expr_df

    #
    # Reorder metadata
    #
    meta_ordered = meta_df.reindex(index=expr_mat.columns)

    batches = meta_ordered[batch_column]

    #
    # ComBat expects:
    # samples x genes
    #

    combat = Combat()
    st.write(batches)

    corrected = combat.fit_transform(expr_mat.T, batches.tolist())

    #
    # Back to:
    # genes x samples
    #
    corrected = corrected.T

    corrected_df = pd.DataFrame(
        corrected, index=expr_mat.index, columns=expr_mat.columns
    )

    return corrected_df


def run_pca(expression_df, apply_log2=False, n_components=2):
    expr = expression_df

    #
    # Force all values to numeric
    # Non-numeric values become NaN
    #
    expr = expr.apply(pd.to_numeric, errors="coerce")

    #
    # Debug information
    #
    total_na = expr.isna().sum().sum()

    print("========== PCA DEBUG ==========")
    print("Expression matrix shape:", expr.shape)
    print("Total missing values:", total_na)

    if total_na > 0:

        print("\nColumns containing NaN:")

        bad_cols = expr.columns[expr.isna().any()]

        print(list(bad_cols))

        print("\nFirst few rows with NaN:")

        print(expr.loc[expr.isna().any(axis=1)].head())

    #
    # Convert to samples x genes
    #
    expr_t = expr.T

    #
    # Optional log2 transform
    #
    if apply_log2:
        expr_t = np.log2(expr_t + 1)

    #
    # Stop immediately if NaN exists
    #
    if expr_t.isna().sum().sum() > 0:

        raise ValueError(
            f"Expression matrix contains "
            f"{expr_t.isna().sum().sum()} "
            f"missing values."
        )

    #
    # Scale
    #
    scaler = StandardScaler()

    scaled_data = scaler.fit_transform(expr_t)

    #
    # PCA
    #
    pca = PCA(n_components=n_components)

    pcs = pca.fit_transform(scaled_data)

    #
    # Output dataframe
    #
    pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_components)])

    pca_df["Sample"] = expr_t.index

    return (pca_df, pca.explained_variance_ratio_)


from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def run_differential_expression(
    expression_df, metadata_df, group_column, group1, group2, apply_log2=False
):
    """
    Differential expression analysis using
    Welch t-test.
    """

    expr = expression_df

    expr = expr.apply(pd.to_numeric, errors="coerce")

    if apply_log2:
        expr = np.log2(expr + 1)

    group1_samples = metadata_df[metadata_df[group_column] == group1].index.tolist()

    group2_samples = metadata_df[metadata_df[group_column] == group2].index.tolist()

    results = []

    for gene in expr.index:

        values1 = expr.loc[gene, group1_samples].dropna()

        values2 = expr.loc[gene, group2_samples].dropna()

        if len(values1) < 2 or len(values2) < 2:
            continue

        mean1 = values1.mean()
        mean2 = values2.mean()

        log2fc = mean1 - mean2

        stat, pvalue = ttest_ind(values1, values2, equal_var=False)

        results.append([gene, mean1, mean2, log2fc, pvalue])
    if len(results) == 0:

        return pd.DataFrame(
            columns=[
                "Gene",
                f"{group1}_Mean",
                f"{group2}_Mean",
                "log2FC",
                "PValue",
                "FDR",
            ]
        )

    de_df = pd.DataFrame(
        results,
        columns=["Gene", f"{group1}_Mean", f"{group2}_Mean", "log2FC", "PValue"],
    )

    de_df["FDR"] = multipletests(de_df["PValue"], method="fdr_bh")[1]

    de_df = de_df.sort_values("FDR")

    return de_df
