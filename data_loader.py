# data_loader.py

import pandas as pd


def load_expression_file(uploaded_file):
    """
    Load expression matrix.

    Expected format:

    Gene,Sample1,Sample2,Sample3
    TP53,10,15,8
    EGFR,20,25,30

    Returns:
        expression_df
    """

    df = pd.read_csv(uploaded_file)

    return df


def validate_expression_matrix(df):
    """
    Validate expression matrix.
    """

    errors = []

    if df.shape[1] < 2:
        errors.append(
            "Expression matrix must contain Gene column and at least one sample column."
        )

    if df.columns[0].lower() != "gene":
        errors.append(
            "First column must be named 'Gene'."
        )

    if df["Gene"].duplicated().any():
        errors.append(
            "Duplicate gene names detected."
        )

    return errors


def summarize_expression_matrix(df):

    n_genes = len(df)

    n_samples = len(df.columns) - 1

    return {
        "Genes": n_genes,
        "Samples": n_samples
    }