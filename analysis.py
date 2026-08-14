# analysis.py

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(expression_df, n_components=2):
    """
    expression_df format:

    Gene  Sample1 Sample2 ...
    TP53  10      15
    EGFR  20      12
    """

    gene_col = expression_df.columns[0]

    expr = expression_df.set_index(gene_col)

    # transpose:
    # samples become rows
    expr_t = expr.T

    scaler = StandardScaler()

    scaled_data = scaler.fit_transform(expr_t)

    pca = PCA(n_components=n_components)

    components = pca.fit_transform(scaled_data)

    result_df = pd.DataFrame(
        components,
        columns=[
            f"PC{i+1}"
            for i in range(n_components)
        ]
    )

    result_df["Sample"] = expr_t.index

    explained_variance = (
        pca.explained_variance_ratio_
    )

    return result_df, explained_variance