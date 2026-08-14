# data_loader.py

import pandas as pd


def load_expression_file(uploaded_file):
    """
    Load expression matrix.

    Expected format:

    Gene,Sample1,Sample2,...
    TP53,10,15,...
    EGFR,20,25,...

    Returns
    -------
    pandas.DataFrame
    """
    return pd.read_csv(uploaded_file)


def load_metadata_file(uploaded_file):
    """
    Load metadata table.

    Expected format:

    Sample,Group,Age,SurvivalTime,Status
    Sample1,Tumor,65,10,1
    Sample2,Normal,70,40,0

    Returns
    -------
    pandas.DataFrame
    """
    return pd.read_csv(uploaded_file)


def validate_expression_matrix(df):
    """
    Validate expression matrix structure.

    Returns
    -------
    list
        List of validation errors
    """

    errors = []

    if df.shape[1] < 2:
        errors.append(
            "Expression matrix must contain a Gene column and at least one sample column."
        )
        return errors

    first_col = df.columns[0]

    if first_col.lower() != "gene":
        errors.append(
            "First column must be named 'Gene'."
        )

    if first_col in df.columns:
        if df[first_col].duplicated().any():
            errors.append(
                "Duplicate gene names detected."
            )

    return errors


def validate_metadata(df):
    """
    Validate metadata structure.

    Returns
    -------
    list
        List of validation errors
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
        