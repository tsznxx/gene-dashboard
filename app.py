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

from analysis import (
    run_pca,
    run_differential_expression
)
from visualization import (
    create_pca_plot,
    create_volcano_plot
)


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

tab_data, tab_pca, tab_de, tab_volcano = st.tabs(
    [
        "Data",
        "PCA",
        "DE Analysis",
        "Volcano Plot"
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
            use_container_width=True,
            config={
                "toImageButtonOptions": {
                    "format": "svg",
                    "filename": "volcano_plot",
                    "height": 800,
                    "width": 1200,
                    "scale": 1
                }
            }
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
        
        
# DE Tab        
with tab_de:

    st.subheader(
        "Differential Expression Analysis"
    )

    eligible_columns = []

    for col in meta_df.columns:

        if meta_df[col].nunique() >= 2:
            eligible_columns.append(col)

    group_column = st.selectbox(
        "Grouping Column",
        eligible_columns,
        key="de_group_column"
    )

    groups = sorted(
        meta_df[group_column]
        .dropna()
        .unique()
        .tolist()
    )

    col1, col2 = st.columns(2)

    with col1:

        group1 = st.selectbox(
            "Group 1",
            groups,
            key="group1"
        )

    with col2:

        remaining = [
            g for g in groups
            if g != group1
        ]

        group2 = st.selectbox(
            "Group 2",
            remaining,
            key="group2"
        )

    if st.button(
        "Run Differential Expression",
        type="primary"
    ):

        de_results = run_differential_expression(
            expression_df=expr_df,
            metadata_df=meta_df,
            group_column=group_column,
            group1=group1,
            group2=group2,
            apply_log2=apply_log2
        )
        st.session_state["de_results"] = de_results

        st.success(
            f"{len(de_results)} genes analysed."
        )

        st.dataframe(
            de_results,
            use_container_width=True
        )

        st.download_button(
            label="Download DE Results",
            data=de_results.to_csv(
                index=False
            ),
            file_name="DE_results.csv",
            mime="text/csv"
        )
        
with tab_volcano:

    st.subheader(
        "Volcano Plot"
    )

    if "de_results" not in st.session_state:

        st.info(
            "Run Differential Expression first."
        )

    else:

        de_df = (
            st.session_state["de_results"]
        )

        significance_column = st.radio(
            "Use significance metric:",
            ["PValue", "FDR"],
            horizontal=True
        )

        col1, col2 = st.columns(2)

        with col1:

            log2fc_cutoff = st.number_input(
                "Absolute log2FC cutoff",
                value=1.0,
                step=0.1
            )

        with col2:

            significance_cutoff = (
                st.number_input(
                    f"{significance_column} cutoff",
                    value=0.05,
                    step=0.01,
                    format="%.3f"
                )
            )
        top_up = (
            de_df
            .sort_values(
                "log2FC",
                ascending=False
            )
            .head(5)["Gene"]
            .tolist()
        )

        top_down = (
            de_df
            .sort_values(
                "log2FC",
                ascending=True
            )
            .head(5)["Gene"]
            .tolist()
        )

        default_genes = (
            top_up + top_down
        )
        custom_gene_text = st.text_area(
            "Additional genes to highlight (comma-separated)",
            value=",".join(default_genes),
            height=100
        )

        highlight_genes = [
            gene.strip()
            for gene in custom_gene_text.split(",")
            if gene.strip()
        ]
        fig = create_volcano_plot(
            de_df=de_df,
            significance_column=
            significance_column,
            significance_cutoff=
            significance_cutoff,
            log2fc_cutoff=
            log2fc_cutoff,
            highlight_genes=
            highlight_genes
        )

        st.plotly_chart(
            fig,
            use_container_width=True,
            config={
                "toImageButtonOptions": {
                    "format": "svg",
                    "filename": "volcano_plot",
                    "height": 800,
                    "width": 1200,
                    "scale": 1
                }
            }
        )
)