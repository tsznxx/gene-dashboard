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
    create_volcano_plot,
    create_gene_boxplot,
    create_heatmap
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

tab_data, tab_pca, tab_de, tab_volcano, tab_box, tab_heatmap = st.tabs(
    [
        "Data",
        "PCA",
        "DE Analysis",
        "Volcano Plot",
        "Gene Boxplot",
        "Heatmap"
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

# ==================================================
# PCA TAB
# ==================================================

with tab_pca:

    st.subheader(
        "Principal Component Analysis"
    )

    color_column = st.selectbox(
        "Color samples by",
        meta_df.columns.tolist(),
        key="pca_color_column"
    )

    st.divider()

    st.subheader(
        "Figure Settings"
    )

    col1, col2 = st.columns(2)

    with col1:

        pca_width = st.number_input(
            "Figure Width (px)",
            min_value=300,
            max_value=1200,
            value=500,
            step=50,
            key="pca_width"
        )

    with col2:

        pca_height = st.number_input(
            "Figure Height (px)",
            min_value=300,
            max_value=1200,
            value=500,
            step=50,
            key="pca_height"
        )

    if st.button(
        "Run PCA",
        key="run_pca_button"
    ):

        try:

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
                pca_df=pca_df,
                color_column=color_column,
                explained_variance=variance,
                width=pca_width,
                height=pca_height
            )

            st.plotly_chart(
                fig,
                width="content",
                config={
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "svg",
                        "filename": "PCA_plot",
                        "width": pca_width,
                        "height": pca_height,
                        "scale": 1
                    }
                }
            )

            st.subheader(
                "Explained Variance"
            )

            st.write(
                {
                    "PC1 (%)": round(
                        variance[0] * 100,
                        2
                    ),
                    "PC2 (%)": round(
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

        except Exception as e:

            st.error(
                str(e)
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

        if de_results.empty:

            st.error(
                "No genes available for differential expression analysis."
            )

            st.stop()
        
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

# ==================================================
# VOLCANO TAB
# ==================================================

with tab_volcano:

    st.subheader("Volcano Plot")

    if "de_results" not in st.session_state:

        st.info(
            "Run Differential Expression first."
        )

    else:

        de_df = st.session_state["de_results"]

        #
        # Volcano Settings
        #
        significance_column = st.radio(
            "Use significance metric",
            ["PValue", "FDR"],
            horizontal=True,
            key="volcano_sig_metric"
        )

        col1, col2 = st.columns(2)

        with col1:
            log2fc_cutoff = st.number_input(
                "Absolute log2FC cutoff",
                value=1.0,
                step=0.1,
                key="volcano_fc_cutoff"
            )

        with col2:
            significance_cutoff = st.number_input(
                f"{significance_column} cutoff",
                value=0.05,
                step=0.01,
                format="%.3f",
                key="volcano_sig_cutoff"
            )

        #
        # Default highlighted genes
        #
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

        default_genes = top_up + top_down

        gene_text = st.text_area(
            "Genes to highlight",
            value=",".join(default_genes),
            height=100,
            key="volcano_gene_text"
        )

        highlight_genes = [
            x.strip()
            for x in gene_text.split(",")
            if x.strip()
        ]

        st.session_state[
            "highlight_genes"
        ] = highlight_genes

        #
        # Figure Settings
        #
        st.divider()
        st.subheader("Figure Settings")

        col3, col4 = st.columns(2)

        with col3:
            volcano_width = st.number_input(
                "Figure Width (px)",
                value=500,
                min_value=200,
                max_value=2000,
                step=100,
                key="volcano_width"
            )

        with col4:
            volcano_height = st.number_input(
                "Figure Height (px)",
                value=500,
                min_value=200,
                max_value=1500,
                step=100,
                key="volcano_height"
            )

        #
        # X-axis
        #
        auto_x = st.checkbox(
            "Automatic X-axis",
            value=True,
            key="volcano_auto_x"
        )

        x_range = None

        if not auto_x:

            c5, c6 = st.columns(2)

            with c5:
                x_min = st.number_input(
                    "X-axis Min",
                    value=-5.0,
                    key="volcano_xmin"
                )

            with c6:
                x_max = st.number_input(
                    "X-axis Max",
                    value=5.0,
                    key="volcano_xmax"
                )

            x_range = [x_min, x_max]

        #
        # Y-axis
        #
        auto_y = st.checkbox(
            "Automatic Y-axis",
            value=True,
            key="volcano_auto_y"
        )

        y_range = None

        if not auto_y:

            c7, c8 = st.columns(2)

            with c7:
                y_min = st.number_input(
                    "Y-axis Min",
                    value=0.0,
                    key="volcano_ymin"
                )

            with c8:
                y_max = st.number_input(
                    "Y-axis Max",
                    value=20.0,
                    key="volcano_ymax"
                )

            y_range = [y_min, y_max]

        #
        # Create Figure
        #
        fig = create_volcano_plot(
            de_df=de_df,
            significance_column=significance_column,
            significance_cutoff=significance_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            highlight_genes=highlight_genes,
            width=volcano_width,
            height=volcano_height,
            x_range=x_range,
            y_range=y_range
        )

        #
        # Display Figure
        #
        st.plotly_chart(
            fig,
            width='content',
            config={
                "displaylogo": False,
                "toImageButtonOptions": {
                    "format": "svg",
                    "filename": "volcano_plot",
                    "width": volcano_width,
                    "height": volcano_height,
                    "scale": 1
                }
            }
        )
        
        
with tab_box:

    st.subheader(
        "Gene Expression Boxplots"
    )

    default_genes = st.session_state.get(
        "highlight_genes",
        []
    )
    st.write(
        st.session_state.get(
            "highlight_genes"
        )
        st.session_state.get(
            "boxplot_gene_text"
        )
        
    )

    gene_text = st.text_area(
        "Genes (comma separated)",
        value=",".join(default_genes),
        height=120,
        key="boxplot_gene_text"
    )

    selected_genes = [

        gene.strip()

        for gene in gene_text.split(",")

        if gene.strip()
    ]

    group_column = st.selectbox(
        "Group By",
        meta_df.columns.tolist(),
        key="boxplot_group_column"
    )

    plot_mode = st.radio(
        "Plot Mode",
        [
            "Combined",
            "Single Gene"
        ],
        horizontal=True
    )

    selected_gene = None

    if (
        plot_mode == "Single Gene"
        and len(selected_genes) > 0
    ):

        selected_gene = st.selectbox(
            "Select Gene",
            selected_genes,
            index=0
        )

    st.divider()

    st.subheader(
        "Figure Settings"
    )

    col1, col2 = st.columns(2)

    with col1:

        plot_width = st.number_input(
            "Figure Width (px)",
            value=1200,
            key="boxplot_width"
        )

    with col2:

        plot_height = st.number_input(
            "Figure Height (px)",
            value=600,
            key="boxplot_height"
        )

    auto_y = st.checkbox(
        "Automatic Y-axis",
        value=True,
        key="boxplot_auto_y"
    )

    y_range = None

    if not auto_y:

        col3, col4 = st.columns(2)

        with col3:

            y_min = st.number_input(
                "Y-axis Minimum",
                value=0.0,
                key="boxplot_ymin"
            )

        with col4:

            y_max = st.number_input(
                "Y-axis Maximum",
                value=20.0,
                key="boxplot_ymax"
            )

        y_range = [
            y_min,
            y_max
        ]

    if st.button(
        "Generate Boxplot",
        key="generate_boxplot"
    ):

        fig = create_gene_boxplot(
            expression_df=expr_df,
            metadata_df=meta_df,
            genes=selected_genes,
            group_column=group_column,
            apply_log2=apply_log2,
            plot_mode=plot_mode,
            selected_gene=selected_gene,
            width=plot_width,
            height=plot_height,
            y_range=y_range
        )

        st.plotly_chart(
            fig,
            width='content',
            config={
                "displaylogo": False,
                "toImageButtonOptions": {
                    "format": "svg",
                    "filename": "gene_boxplot",
                    "width": plot_width,
                    "height": plot_height,
                    "scale": 1
                }
            }
            
        )
        
# ==================================================
# HEATMAP TAB
# ==================================================

with tab_heatmap:

    st.subheader(
        "Expression Heatmap"
    )

    default_genes = st.session_state.get(
        "highlight_genes",
        []
    )

    gene_text = st.text_area(
        "Genes (comma separated)",
        value=",".join(default_genes),
        height=120,
        key="heatmap_genes"
    )

    selected_genes = [

        x.strip()

        for x in gene_text.split(",")

        if x.strip()
    ]

    annotation_column = (
        st.selectbox(
            "Annotation Column",
            meta_df.columns,
            index=1,
            key="heatmap_annotation"
        )
    )

    col1, col2 = st.columns(2)

    with col1:

        zscore_by_gene = (
            st.checkbox(
                "Z-score by Gene",
                value=True,
                key="heatmap_zscore"
            )
        )

        cluster_genes = (
            st.checkbox(
                "Cluster Genes",
                value=True,
                key="heatmap_cluster_genes"
            )
        )

    with col2:

        cluster_samples = (
            st.checkbox(
                "Cluster Samples",
                value=True,
                key="heatmap_cluster_samples"
            )
        )

        colorscale = (
            st.selectbox(
                "Color Map",
                [
                    "RdBu_r",
                    "Viridis",
                    "Plasma",
                    "Magma"
                ],
                key="heatmap_colorscale"
            )
        )

    st.divider()

    st.subheader(
        "Figure Settings"
    )

    c1, c2 = st.columns(2)

    with c1:

        heatmap_width = st.number_input(
            "Figure Width (px)",
            value=1200,
            min_value=300,
            max_value=4000,
            step=100,
            key="heatmap_width"
        )

    with c2:

        heatmap_height = (
            st.number_input(
                "Figure Height (px)",
                value=800,
                min_value=300,
                max_value=4000,
                step=100,
                key="heatmap_height"
            )
        )
    if len(selected_genes) == 0:

        st.info(
            "Run DE analysis first or provide at least one gene."
        )

    else:
        fig = create_heatmap(
            expression_df=expr_df,
            metadata_df=meta_df,
            genes=selected_genes,
            annotation_column=
            annotation_column,
            apply_log2=apply_log2,
            zscore_by_gene=
            zscore_by_gene,
            cluster_samples=
            cluster_samples,
            cluster_genes=
            cluster_genes,
            width=heatmap_width,
            height=heatmap_height,
            colorscale=colorscale
        )

        st.plotly_chart(
            fig,
            width="content",
            config={
                "displaylogo": False,
                "toImageButtonOptions": {
                    "format": "svg",
                    "filename":"expression_heatmap",
                    "width":heatmap_width,
                    "height": plot_height,
                    "scale": 1
                }
            }       
        )