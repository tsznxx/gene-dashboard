import streamlit as st

from data_loader import (
    load_expression_file,
    load_metadata_file,
    validate_expression_matrix,
    validate_metadata,
    validate_sample_matching,
    summarize_expression,
    summarize_metadata
)

from analysis import run_pca
from visualization import create_pca_plot


st.set_page_config(
    page_title="Gene Expression Dashboard",
    layout="wide"
)

st.title("Gene Expression Dashboard")

#
# Sidebar
#

with st.sidebar:

    st.header("Data Upload")

    uploaded_expression = st.file_uploader(
        "Expression Matrix (CSV)",
        type=["csv"]
    )

    uploaded_metadata = st.file_uploader(
        "Metadata Table (CSV)",
        type=["csv"]
    )

#
# Main
#

if not uploaded_expression or not uploaded_metadata:

    st.info(
        "Upload both expression matrix and metadata table to begin."
    )

    st.stop()

#
# load data
#

expr_df = load_expression_file(
    uploaded_expression
)

meta_df = load_metadata_file(
    uploaded_metadata
)

expr_errors = validate_expression_matrix(
    expr_df
)

meta_errors = validate_metadata(
    meta_df
)

all_errors = expr_errors + meta_errors

if all_errors:

    st.error(
        "Validation failed."
    )

    for err in all_errors:
        st.write(f"• {err}")

    st.stop()

#
# sample matching
#

matching_result = validate_sample_matching(
    expr_df,
    meta_df
)

if not matching_result["matching"]:

    st.error(
        "Sample mismatch detected."
    )

    if matching_result["missing_in_metadata"]:

        st.write(
            "Missing in metadata:"
        )

        st.write(
            matching_result["missing_in_metadata"]
        )

    if matching_result["missing_in_expression"]:

        st.write(
            "Missing in expression matrix:"
        )

        st.write(
            matching_result["missing_in_expression"]
        )

    st.stop()

st.success(
    "Data loaded successfully."
)

#
# save to session state
#

st.session_state["expression_df"] = expr_df
st.session_state["metadata_df"] = meta_df

#
# Global options
#

with st.sidebar:

    st.header("Global Settings")

    apply_log2 = st.checkbox(
        "Apply log2(x+1) transformation",
        value=True,
        help="Recommended for RNA-seq count data"
    )

#
# Tabs
#

tab_data, tab_pca = st.tabs(
    [
        "Data",
        "PCA"
    ]
)

#
# Data Tab
#

with tab_data:

    st.subheader("Dataset Summary")

    expr_summary = summarize_expression(
        expr_df
    )

    meta_summary = summarize_metadata(
        meta_df
    )

    c1, c2, c3 = st.columns(3)

    with c1:
        st.metric(
            "Genes",
            expr_summary["Genes"]
        )

    with c2:
        st.metric(
            "Samples",
            expr_summary["Samples"]
        )

    with c3:
        st.metric(
            "Groups",
            meta_summary.get(
                "Groups",
                "NA"
            )
        )

    st.divider()

    st.subheader(
        "Expression Matrix Preview"
    )

    st.dataframe(
        expr_df.head(20),
        use_container_width=True
    )

    st.subheader(
        "Metadata Preview"
    )

    st.dataframe(
        meta_df.head(20),
        use_container_width=True
    )

#
# PCA Tab
#

with tab_pca:

    st.subheader(
        "Principal Component Analysis"
    )

    color_column = st.selectbox(
        "Color samples by",
        meta_df.columns.tolist()
    )

    if st.button(
        "Run PCA",
        type="primary"
    ):

        pca_df, variance = run_pca(
            expr_df,
            apply_log2=apply_log2
        )

        pca_df = pca_df.merge(
            meta_df,
            on="Sample",
            how="left"
        )

        fig = create_pca_plot(
            pca_df,
            color_column,
            variance
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

        st.subheader(
            "Explained Variance"
        )

        st.write(
            {
                "PC1": round(
                    variance[0] * 100,
                    2
                ),
                "PC2": round(
                    variance[1] * 100,
                    2
                )
            }
        )

        st.subheader(
            "PCA Coordinates"
        )

        st.dataframe(
            pca_df,
            use_container_width=True
        )