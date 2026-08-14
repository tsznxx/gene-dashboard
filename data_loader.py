# data_loader.py

import pandas as pd
import numpy as np


def load_expression_file(uploaded_file):
    """
    Load expression matrix.

    Expected format:

    Gene,Sample1,Sample2,...
    TP53,10,15,...
    EGFR,20,25,...
    """

    return pd.read_csv(uploaded_file)


def load_metadata_file(uploaded_file):
    """
    Load metadata file.

    Expected format:

    Sample,Group,Age,SurvivalTime,Status
    Sample1,Tumor,65,10,1
    Sample2,Normal,70,40,0
    """

    return pd.read_csv(uploaded_file)


def validate_expression_matrix(df):
    """
    Validate expression matrix.
    """

    errors = []

    if df.shape[1] < 2:
        errors.append(
            "Expression matrix must contain Gene column and at least one sample column."
        )
        return errors

    gene_col = df.columns[0]

    if gene_col.lower() != "gene":
        errors.append(
            "First column must be named 'Gene'."
        )

    if df[gene_col].duplicated().any():
        errors.append(
            "Duplicate gene names detected."
        )

    return errors


def validate_metadata(df):
    """
    Validate metadata table.
    """

    errors = []

    if "Sample" not in df.columns:
        errors.append(
            "Metadata must contain a 'Sample' column."
        )
        return errors

    if df["Sample"].duplicated().any():
        errors.append(
            "Duplicate Sample IDs detected."
        )

    return errors


def validate_sample_matching(expr_df, meta_df):
    """
    Check if sample IDs match between
    expression matrix and metadata.
    """

    expr_samples = set(expr_df.columns[1:])

    meta_samples = set(meta_df["Sample"])

    missing_in_metadata = expr_samples - meta_samples
    missing_in_expression = meta_samples - expr_samples

    return {
        "matching": (
            len(missing_in_metadata) == 0
            and len(missing_in_expression) == 0
        ),
        "missing_in_metadata": sorted(
            list(missing_in_metadata)
        ),
        "missing_in_expression": sorted(
            list(missing_in_expression)
        )
    }


def summarize_expression(expr_df):
    """
    Summarize expression matrix.
    """

    return {
        "Genes": expr_df.shape[0],
        "Samples": expr_df.shape[1] - 1
    }


def summarize_metadata(meta_df):
    """
    Summarize metadata.
    """

    summary = {
        "Samples": len(meta_df)
    }

    if "Group" in meta_df.columns:
        summary["Groups"] = (
            meta_df["Group"].nunique()
        )

    return summary


def get_gene_names(expr_df):
    """
    Return gene names.
    """

    return expr_df.iloc[:, 0].tolist()


def get_sample_names(expr_df):
    """
    Return sample names.
    """

    return expr_df.columns[1:].tolist()


def prepare_expression_matrix(
    expr_df,
    apply_log2=False
):
    """
    Convert expression matrix into a
    samples x genes numeric matrix.

    Output format:

               TP53  EGFR  MYC
    Sample1     10    20     5
    Sample2     15    22     3
    """

    gene_col = expr_df.columns[0]

    matrix = expr_df.set_index(
        gene_col
    )

    matrix = matrix.astype(float)

    if apply_log2:
        matrix = np.log2(matrix + 1)

    # transpose:
    # rows=samples, columns=genes
    matrix = matrix.T

    return matrix