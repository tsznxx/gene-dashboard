# app.py

import streamlit as st
import pandas as pd

from data_loader import (
    load_expression_file,
    validate_expression_matrix,
    summarize_expression_matrix
)

st.set_page_config(
    page_title="Gene Expression Dashboard",
    layout="wide"
)

st.title("Gene Expression Dashboard")

st.markdown(
    """
    ### Module 1: Upload Expression Matrix

    Upload a gene expression matrix where:

    - Rows = Genes
    - Columns = Samples
    - First column = Gene
    """
)

uploaded_expression = st.file_uploader(
    "Upload Expression Matrix (.csv)",
    type=["csv"]
)

if uploaded_expression:

    try:

        expression_df = load_expression_file(
            uploaded_expression
        )

        errors = validate_expression_matrix(
            expression_df
        )

        if errors:

            st.error("Validation failed")

            for err in errors:
                st.write("•", err)

        else:

            st.success("Expression matrix loaded successfully")

            summary = summarize_expression_matrix(
                expression_df
            )

            col1, col2 = st.columns(2)

            with col1:
                st.metric("Genes", summary["Genes"])

            with col2:
                st.metric("Samples", summary["Samples"])

            st.subheader("Preview")

            st.dataframe(
                expression_df.head(20),
                use_container_width=True
            )

    except Exception as e:

        st.exception(e)