# analysis.py

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(
    expression_df,
    apply_log2=False,
    n_components=2
):

    gene_col = expression_df.columns[0]

    #
    # Set Gene as row index
    #
    expr = expression_df.set_index(gene_col)

    #
    # Force all values to numeric
    # Non-numeric values become NaN
    #
    expr = expr.apply(
        pd.to_numeric,
        errors="coerce"
    )

    #
    # Debug information
    #
    total_na = expr.isna().sum().sum()

    print("========== PCA DEBUG ==========")
    print("Expression matrix shape:", expr.shape)
    print("Total missing values:", total_na)

    if total_na > 0:

        print("\nColumns containing NaN:")

        bad_cols = expr.columns[
            expr.isna().any()
        ]

        print(list(bad_cols))

        print("\nFirst few rows with NaN:")

        print(
            expr.loc[
                expr.isna().any(axis=1)
            ].head()
        )

    #
    # Convert to samples x genes
    #
    expr_t = expr.T

    #
    # Optional log2 transform
    #
    if apply_log2:
        expr_t = np.log2(
            expr_t + 1
        )

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

    scaled_data = scaler.fit_transform(
        expr_t
    )

    #
    # PCA
    #
    pca = PCA(
        n_components=n_components
    )

    pcs = pca.fit_transform(
        scaled_data
    )

    #
    # Output dataframe
    #
    pca_df = pd.DataFrame(
        pcs,
        columns=[
            f"PC{i+1}"
            for i in range(
                n_components
            )
        ]
    )

    pca_df["Sample"] = expr_t.index

    return (
        pca_df,
        pca.explained_variance_ratio_
    )