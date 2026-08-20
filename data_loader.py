# data_loader.py

import pandas as pd
import numpy as np


def load_table(file):

    if file.name.endswith(".csv"):
        df = pd.read_csv(file,index_col=0)
    elif (
        file.name.endswith(".tsv")
        or file.name.endswith(".txt")
    ):

        df =  pd.read_csv(file,sep="\t",index_col=0)

    else:

        raise ValueError(
            "Unsupported file type."
        )
    df.index = df.index.str.strip()
    df.columns = df.columns.str.strip()    

def load_expression_file(uploaded_file):
    """
    Load expression matrix.

    Expected format:

    Gene,Sample1,Sample2,...
    TP53,10,15,...
    EGFR,20,25,...
    """

    return load_table(uploaded_file)


def load_metadata_file(uploaded_file):
    """
    Load metadata file.

    Expected format:

    Sample,Group,Age,SurvivalTime,Status
    Sample1,Tumor,65,10,1
    Sample2,Normal,70,40,0
    """

    return load_table(uploaded_file)


def validate_expression_matrix(df):
    """
    Validate expression matrix.
    """

    errors = []

    if df.shape[1] < 4:
        errors.append(
            "Expression matrix must contain Gene column and at least three sample columns."
        )
        return errors

    if df.index.duplicated().any():
        errors.append(
            "Duplicate gene names detected."
        )

    return errors


def validate_metadata(df):
    """
    Validate metadata table.
    """

    errors = []

    if df.index.duplicated().any():
        errors.append(
            "Duplicate Sample IDs detected."
        )

    return errors


def validate_sample_matching(expr_df, meta_df):
    """
    Check if sample IDs match between
    expression matrix and metadata.
    """

    expr_samples = set(expr_df.columns.to_list())

    meta_samples = set(meta_df.index.to_list())

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
        "Samples": expr_df.shape[1]
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

    return expr_df.index.tolist()


def get_sample_names(expr_df):
    """
    Return sample names.
    """

    return expr_df.columns.tolist()


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

    matrix = expr_df

    matrix = matrix.astype(float)

    if apply_log2:
        matrix = np.log2(matrix + 1)

    # transpose:
    # rows=samples, columns=genes
    matrix = matrix.T

    return matrix