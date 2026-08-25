# analysis.py

# for test
import streamlit as st

import numpy as np
import pandas as pd
from pycombat import Combat
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

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

    corrected = combat.fit_transform(expr_mat.T.values, batches.values)

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





def run_differential_expression(
    expression_df, metadata_df, group_column, group1, group2, apply_log2=False
):
    """
    Differential expression analysis using
    Welch t-test.
    """

    expr = expression_df

    expr = expr.apply(pd.to_numeric, errors="coerce")

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


# ==================================================
# Gene Validation
# ==================================================

def validate_correlation_genes(
    expr_df,
    genes
):
    """
    Validate a gene list against the expression matrix.

    Parameters
    ----------
    expr_df : pandas.DataFrame
        Expression matrix with genes as rows and samples
        as columns.

    genes : list
        User-provided gene list.

    Returns
    -------
    valid_genes : list
        Genes found in the expression matrix.

    missing_genes : list
        Genes not found in the expression matrix.
    """

    if genes is None:
        genes = []

    cleaned_genes = []

    for gene in genes:

        gene = str(gene).strip()

        if (
            gene
            and gene not in cleaned_genes
        ):
            cleaned_genes.append(gene)

    valid_genes = [
        gene
        for gene in cleaned_genes
        if gene in expr_df.index
    ]

    missing_genes = [
        gene
        for gene in cleaned_genes
        if gene not in expr_df.index
    ]

    return valid_genes, missing_genes


# ==================================================
# Signature Calculation
# ==================================================

def calculate_signature_score(
    expr_df,
    genes,
    aggregation="Mean"
):
    """
    Calculate one sample-level score from a gene signature.

    Unlike a Gene List, a Gene Signature is aggregated
    into one vector.

    Parameters
    ----------
    expr_df : pandas.DataFrame
        Genes x samples expression matrix.

    genes : list
        Genes included in the signature.

    aggregation : str
        Mean, Median, Sum, or Mean Z-score.

    Returns
    -------
    pandas.Series
        Signature score for each sample.
    """

    valid_genes, missing_genes = (
        validate_correlation_genes(
            expr_df,
            genes
        )
    )

    if len(valid_genes) == 0:

        raise ValueError(
            "None of the signature genes were found "
            "in the expression matrix."
        )

    signature_matrix = (
        expr_df.loc[
            valid_genes
        ]
        .apply(
            pd.to_numeric,
            errors="coerce"
        )
    )

    if aggregation == "Mean":

        score = signature_matrix.mean(
            axis=0
        )

    elif aggregation == "Median":

        score = signature_matrix.median(
            axis=0
        )

    elif aggregation == "Sum":

        score = signature_matrix.sum(
            axis=0
        )

    elif aggregation == "Mean Z-score":

        row_mean = signature_matrix.mean(
            axis=1
        )

        row_std = signature_matrix.std(
            axis=1,
            ddof=1
        ).replace(
            0,
            np.nan
        )

        zscore_matrix = (
            signature_matrix
            .sub(
                row_mean,
                axis=0
            )
            .div(
                row_std,
                axis=0
            )
        )

        score = zscore_matrix.mean(
            axis=0
        )

    else:

        raise ValueError(
            f"Unsupported signature aggregation: "
            f"{aggregation}"
        )

    score = pd.to_numeric(
        score,
        errors="coerce"
    )

    return score


# ==================================================
# Subject Resolution
# ==================================================

def resolve_correlation_subject(
    expr_df,
    subject_type,
    gene=None,
    gene_list=None,
    signature_name=None,
    signature_aggregation="Mean"
):
    """
    Resolve a correlation subject into named vectors.

    Every returned item is a separate sample-level vector.

    Subject behavior
    ----------------
    Single Gene:
        One named vector.

    Gene List:
        One separate vector per gene. The genes are not
        aggregated.

    Gene Signature:
        One aggregated vector with a user-provided name.

    All Genes:
        One separate vector per gene in the expression
        matrix.

    Returns
    -------
    dict
        {
            "subject_type": str,
            "vectors": {
                subject_name: pandas.Series,
                ...
            },
            "valid_genes": list,
            "missing_genes": list
        }
    """

    supported_types = [
        "Single Gene",
        "Gene List",
        "Gene Signature",
        "All Genes"
    ]

    if subject_type not in supported_types:

        raise ValueError(
            f"Unsupported subject type: {subject_type}"
        )

    # ----------------------------------
    # Single Gene
    # ----------------------------------

    if subject_type == "Single Gene":

        if gene is None:

            raise ValueError(
                "A gene must be selected."
            )

        if gene not in expr_df.index:

            raise ValueError(
                f"Gene '{gene}' was not found "
                "in the expression matrix."
            )

        vector = pd.to_numeric(
            expr_df.loc[gene],
            errors="coerce"
        )

        return {
            "subject_type": subject_type,
            "vectors": {
                str(gene): vector
            },
            "valid_genes": [
                str(gene)
            ],
            "missing_genes": []
        }

    # ----------------------------------
    # Gene List
    # ----------------------------------

    if subject_type == "Gene List":

        valid_genes, missing_genes = (
            validate_correlation_genes(
                expr_df,
                gene_list
            )
        )

        if len(valid_genes) == 0:

            raise ValueError(
                "None of the genes in the gene list "
                "were found in the expression matrix."
            )

        vectors = {}

        for current_gene in valid_genes:

            vectors[
                current_gene
            ] = pd.to_numeric(
                expr_df.loc[
                    current_gene
                ],
                errors="coerce"
            )

        return {
            "subject_type": subject_type,
            "vectors": vectors,
            "valid_genes": valid_genes,
            "missing_genes": missing_genes
        }

    # ----------------------------------
    # Gene Signature
    # ----------------------------------

    if subject_type == "Gene Signature":

        valid_genes, missing_genes = (
            validate_correlation_genes(
                expr_df,
                gene_list
            )
        )

        if len(valid_genes) == 0:

            raise ValueError(
                "None of the signature genes were found "
                "in the expression matrix."
            )

        clean_signature_name = (
            str(signature_name).strip()
            if signature_name is not None
            else ""
        )

        if not clean_signature_name:

            clean_signature_name = (
                "Gene Signature"
            )

        signature_score = (
            calculate_signature_score(
                expr_df=expr_df,
                genes=valid_genes,
                aggregation=
                signature_aggregation
            )
        )

        return {
            "subject_type": subject_type,
            "vectors": {
                clean_signature_name:
                signature_score
            },
            "valid_genes": valid_genes,
            "missing_genes": missing_genes
        }

    # ----------------------------------
    # All Genes
    # ----------------------------------

    numeric_expr = expr_df.apply(
        pd.to_numeric,
        errors="coerce"
    )

    vectors = {
        str(current_gene):
        numeric_expr.loc[current_gene]
        for current_gene
        in numeric_expr.index
    }

    return {
        "subject_type": subject_type,
        "vectors": vectors,
        "valid_genes": (
            numeric_expr.index
            .astype(str)
            .tolist()
        ),
        "missing_genes": []
    }


# ==================================================
# Pairwise Correlation
# ==================================================

def compute_pairwise_correlation(
    x,
    y,
    method="pearson",
    min_samples=3
):
    """
    Calculate correlation between two sample-level vectors.

    Missing values are removed pairwise.

    Parameters
    ----------
    x, y : pandas.Series
        Sample-level vectors.

    method : str
        pearson or spearman.

    min_samples : int
        Minimum number of paired observations.

    Returns
    -------
    coefficient : float
    pvalue : float
    n_samples : int
    """

    pair_df = pd.concat(
        [
            pd.to_numeric(
                x,
                errors="coerce"
            ).rename("x"),
            pd.to_numeric(
                y,
                errors="coerce"
            ).rename("y")
        ],
        axis=1,
        join="inner"
    )

    pair_df = (
        pair_df
        .replace(
            [np.inf, -np.inf],
            np.nan
        )
        .dropna()
    )

    n_samples = len(
        pair_df
    )

    if n_samples < min_samples:

        return (
            np.nan,
            np.nan,
            n_samples
        )

    # Correlation is undefined for a constant vector.
    if (
        pair_df["x"].nunique() < 2
        or pair_df["y"].nunique() < 2
    ):

        return (
            np.nan,
            np.nan,
            n_samples
        )

    method = method.lower()

    if method == "pearson":

        coefficient, pvalue = pearsonr(
            pair_df["x"],
            pair_df["y"]
        )

    elif method == "spearman":

        coefficient, pvalue = spearmanr(
            pair_df["x"],
            pair_df["y"]
        )

    else:

        raise ValueError(
            f"Unsupported correlation method: "
            f"{method}"
        )

    return (
        float(coefficient),
        float(pvalue),
        n_samples
    )


# ==================================================
# FDR Correction
# ==================================================

def add_correlation_fdr(
    corr_df
):
    """
    Add Benjamini-Hochberg FDR values.

    Invalid or undefined p-values remain NaN.
    """

    result_df = corr_df.copy()

    result_df["FDR"] = np.nan

    valid_mask = (
        result_df["PValue"]
        .notna()
        &
        np.isfinite(
            result_df["PValue"]
        )
    )

    if valid_mask.any():

        corrected_values = (
            multipletests(
                result_df.loc[
                    valid_mask,
                    "PValue"
                ],
                method="fdr_bh"
            )[1]
        )

        result_df.loc[
            valid_mask,
            "FDR"
        ] = corrected_values

    return result_df


# ==================================================
# Main Correlation Analysis
# ==================================================

def run_correlation_analysis(
    subject_a,
    subject_b,
    method="pearson",
    min_samples=3,
    skip_identical_pairs=True,
    allow_all_vs_all=False
):
    """
    Perform pairwise correlation across two resolved subjects.

    Examples
    --------
    Single Gene vs Single Gene:
        1 correlation.

    Single Gene vs Gene List:
        One correlation per gene in the list.

    Gene List vs Gene List:
        All pairwise correlations between the two lists.

    Gene Signature vs All Genes:
        Signature score against each gene.

    All Genes vs All Genes:
        Disabled by default because it can generate a very
        large number of comparisons.

    Parameters
    ----------
    subject_a, subject_b : dict
        Output from resolve_correlation_subject().

    method : str
        pearson or spearman.

    min_samples : int
        Minimum number of paired observations.

    skip_identical_pairs : bool
        Skip a gene correlated with itself when both subjects
        contain the same named gene.

    allow_all_vs_all : bool
        Whether to allow All Genes vs All Genes.

    Returns
    -------
    pandas.DataFrame
        Subject_A, Subject_B, Coefficient, PValue, FDR,
        and N_Samples.
    """

    type_a = subject_a[
        "subject_type"
    ]

    type_b = subject_b[
        "subject_type"
    ]

    if (
        type_a == "All Genes"
        and type_b == "All Genes"
        and not allow_all_vs_all
    ):

        raise ValueError(
            "All Genes vs All Genes is disabled because "
            "the number of pairwise comparisons can be "
            "extremely large. Please use a gene list or "
            "gene signature for at least one subject."
        )

    vectors_a = subject_a[
        "vectors"
    ]

    vectors_b = subject_b[
        "vectors"
    ]

    if len(vectors_a) == 0:

        raise ValueError(
            "Subject A contains no valid vectors."
        )

    if len(vectors_b) == 0:

        raise ValueError(
            "Subject B contains no valid vectors."
        )

    results = []

    for name_a, vector_a in (
        vectors_a.items()
    ):

        for name_b, vector_b in (
            vectors_b.items()
        ):

            if (
                skip_identical_pairs
                and name_a == name_b
            ):
                continue

            coefficient, pvalue, n_samples = (
                compute_pairwise_correlation(
                    x=vector_a,
                    y=vector_b,
                    method=method,
                    min_samples=min_samples
                )
            )

            results.append(
                {
                    "Subject_A":
                    name_a,

                    "Subject_B":
                    name_b,

                    "Coefficient":
                    coefficient,

                    "PValue":
                    pvalue,

                    "N_Samples":
                    n_samples
                }
            )

    corr_df = pd.DataFrame(
        results,
        columns=[
            "Subject_A",
            "Subject_B",
            "Coefficient",
            "PValue",
            "N_Samples"
        ]
    )

    if corr_df.empty:

        corr_df["FDR"] = pd.Series(
            dtype=float
        )

        return corr_df

    corr_df = add_correlation_fdr(
        corr_df
    )

    corr_df["Absolute_Correlation"] = (
        corr_df["Coefficient"].abs()
    )

    corr_df = (
        corr_df
        .sort_values(
            by=[
                "Absolute_Correlation",
                "PValue"
            ],
            ascending=[
                False,
                True
            ],
            na_position="last"
        )
        .reset_index(
            drop=True
        )
    )

    return corr_df


# ==================================================
# Retrieve Vectors for Scatter Plot
# ==================================================

def get_correlation_plot_vectors(
    subject_a,
    subject_b,
    subject_a_name,
    subject_b_name
):
    """
    Retrieve the exact vectors for a selected correlation
    result row.

    This is useful for creating the correlation scatter plot.
    """

    vectors_a = subject_a[
        "vectors"
    ]

    vectors_b = subject_b[
        "vectors"
    ]

    if subject_a_name not in vectors_a:

        raise ValueError(
            f"Subject A '{subject_a_name}' "
            "is unavailable."
        )

    if subject_b_name not in vectors_b:

        raise ValueError(
            f"Subject B '{subject_b_name}' "
            "is unavailable."
        )

    x = vectors_a[
        subject_a_name
    ]

    y = vectors_b[
        subject_b_name
    ]

    plot_df = pd.concat(
        [
            pd.to_numeric(
                x,
                errors="coerce"
            ).rename(
                subject_a_name
            ),
            pd.to_numeric(
                y,
                errors="coerce"
            ).rename(
                subject_b_name
            )
        ],
        axis=1,
        join="inner"
    )

    plot_df = (
        plot_df
        .replace(
            [np.inf, -np.inf],
            np.nan
        )
        .dropna()
    )

    plot_df.index.name = "Sample"

    return plot_df.reset_index()