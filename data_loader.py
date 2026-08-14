# data_loader.py

import pandas as pd


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
    Load metadata table.

    Expected format:

    Sample,Group,Age,SurvivalTime,Status
    Sample1,Tumor,65,10,1
    Sample2,Normal,70,40,0
    """

    return pd.read_csv(uploaded_file)


def validate_expression_matrix(df):
    """
    Validate expression matrix structure.
    """

    errors = []

    if df.shape[1] < 2:
        errors.append(
            "Expression matrix must contain a Gene column and at least one sample column."
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
    Validate metadata structure.
    """

    errors = []

    sample_col = None

    for col in df.columns:
        if col.lower() == "sample":
            sample_col = col
            break

    if sample_col is None:
        errors.append(
            "Metadata table must contain a 'Sample' column."
        )
        return errors

    if df[sample_col].duplicated().any():
        errors.append(
            "Duplicate sample IDs detected."
        )

    return errors


def validate_sample_matching(expr_df, meta_df):
    """
    Check whether samples match between
    expression matrix and metadata table.
    """

    expr_samples = set(expr_df.columns[1:])

    sample_col = [
        col
        for col in meta_df.columns
        if col.lower() == "sample"
    ][0]

    meta_samples = set(meta_df[sample_col])

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
    Generate summary statistics
    for expression matrix.
    """

    return {
        "Genes": expr_df.shape[0],
        "Samples": expr_df.shape[1] - 1
    }


def summarize_metadata(meta_df):
    """
    Generate summary statistics
    for metadata table.
    """

    summary = {
        "Samples": len(meta_df)
    }

    group_cols = [
        col
        for col in meta_df.columns
        if col.lower() == "group"
    ]

    if len(group_cols) > 0:
        summary["Groups"] = (
            meta_df[group_cols[0]]
            .nunique()
        )

    return summary


def get_sample_column(meta_df):
    """
    Return metadata sample column name.
    """

    for col in meta_df.columns:
        if col.lower() == "sample"