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
    """
    Run PCA on expression matrix.

    Parameters
    ----------
    expression_df : DataFrame
        Expression matrix:
        Gene, Sample1, Sample2 ...

    apply_log2 : bool
        Apply log2(x+1) transformation

    n_components : int
        Number of principal components

    Returns
    -------
    pca_df : DataFrame
    explained_variance : ndarray
    """

    gene_col = expression_df.columns[0]

    expr = expression_df.set_index(gene_col)

    # genes x samples -> samples x genes
    expr_t = expr.T.astype(float)

    if apply_log2:
        expr_t = np.log2(expr_t + 1)

    scaler = StandardScaler()

    scaled_data = scaler.fit_transform(expr_t)

    pca = PCA(n_components=n_components)

    pcs = pca.fit_transform(scaled_data)

    pca_df = pd.DataFrame(
        pcs,
        columns=[
            f"PC{i+1}"
            for i in range(n_components)
        ]
    )

    pca_df["Sample"] = expr_t.index

    return (
        pca_df,
        pca.explained_variance_ratio_
    )