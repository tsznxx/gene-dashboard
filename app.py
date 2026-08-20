import streamlit as st
import pandas as pd

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

st.title("Gene Expression Dashboard")

st.set_page_config(
    page_title="Gene Expression Dashboard",
    page_icon="assets/logo.png",
    layout="wide"
)

st.markdown(
    """
    <style>
    .footer {
        position: fixed;
        left: 0;
        bottom: 0;
        width: 100%;
        background-color: white;
        color: gray;
        text-align: center;
        font-size: 12px;
        padding: 5px;
        border-top: 1px solid #dddddd;
        z-index: 999;
    }
    </style>

    <div class="footer">
        © 2026 H. Lee Moffitt Cancer Center | Gene Expression Dashboard
    </div>
    """,
    unsafe_allow_html=True
)

#
# Sidebar
#

# ==================================================
# SIDEBAR
# ==================================================

with st.sidebar:

    st.image(
        "assets/logo.png",
        width='content'
    )

    st.divider()

    st.header("Data Source")

    #
    # Load Example Dataset
    #
    if st.button(
        "Load Example Dataset",
        type="primary",
        key="load_example_dataset"
    ):

        expr_df_example = pd.read_csv(
            "example_data/example_expression.tsv",sep='\t',index_col=0
        )

        meta_df_example = pd.read_csv(
            "example_data/example_metadata.tsv",sep='\t',index_col=0
        )

        st.session_state[
            "expr_df"
        ] = expr_df_example

        st.session_state[
            "meta_df"
        ] = meta_df_example

        st.session_state[
            "using_example_data"
        ] = True


    if st.session_state.get(
        "using_example_data",
        False
    ):

        st.success(
            "Using example dataset"
        )

    st.divider()

    #
    # Upload Data
    #
    uploaded_expression = st.file_uploader(
        "Expression Matrix",
        type=[
            "csv",
            "tsv",
            "txt"
        ],
        key="expression_upload"
    )

    uploaded_metadata = st.file_uploader(
        "Metadata Table",
        type=[
            "csv",
            "tsv",
            "txt"
        ],
        key="metadata_upload"
    )

    #
    # Global Settings
    #
    st.divider()

    st.header(
        "Global Settings"
    )

    apply_log2 = st.checkbox(
        "Apply log2(x+1)",
        value=True,
        key="apply_log2"
    )

    #
    # Example Files
    #
    st.divider()

    st.header(
        "Example Files"
    )

    with open(
        "example_data/example_expression.tsv",
        "rb"
    ) as f:

        st.download_button(
            "Download Example Expression",
            data=f,
            file_name=
            "example_expression.tsv",
            mime="text/tab-separated-values",
            key="download_example_expr"
        )

    with open(
        "example_data/example_metadata.tsv",
        "rb"
    ) as f:

        st.download_button(
            "Download Example Metadata",
            data=f,
            file_name=
            "example_metadata.tsv",
            mime="text/tab-separated-values",
            key="download_example_meta"
        )

        
#
# Main
#


# ==================================================
# LOAD DATA
# ==================================================

expr_df = None
meta_df = None

#
# Example Dataset
#
if st.session_state.get(
    "using_example_data",
    False
):

    expr_df = st.session_state[
        "expr_df"
    ]

    meta_df = st.session_state[
        "meta_df"
    ]

#
# Uploaded Dataset
#
elif (
    uploaded_expression is not None
    and uploaded_metadata is not None
):

    expr_df = load_expression_file(
        uploaded_expression
    )

    meta_df = load_metadata_file(
        uploaded_metadata
    )

#
# Nothing Loaded Yet
#
else:

    st.info(
        "Upload files or click "
        "'Load Example Dataset'."
    )

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
        width='content'
    )

    st.subheader(
        "Metadata Preview"
    )

    st.dataframe(
        meta_df.head(20),
        width='content'
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
        meta_df.columns,
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
            value=450,
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
                width='content'
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
        de_results = de_results.sort_values('log2FC')
        st.session_state["de_results"] = de_results
        all_genes = sorted(de_results["Gene"].astype(str).unique().tolist())
        st.session_state["all_genes"] = all_genes

        st.success(
            f"{len(de_results)} genes analysed."
        )

        st.dataframe(
            de_results,
            width='content'
        )

        st.download_button(
            label="Download DE Results",
            data=de_results.to_csv(
                index=False,sep='\t'
            ),
            file_name="DE_results.tsv",
            mime="text/tab-separated-values"
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

        # --------------------------------------------------
        # Volcano Settings
        # --------------------------------------------------

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
                key="volcano_sig_cutoff"
            )

        # --------------------------------------------------
        # Highlighted Genes
        # --------------------------------------------------

        st.subheader(
            "Genes to Highlight"
        )

        use_top_genes = st.checkbox(
            "Use Top Significant Genes",
            value=True,
            key="volcano_use_top_genes"
        )

        if "highlight_genes" not in st.session_state:

            st.session_state["highlight_genes"] = []

        if "volcano_gene_text" not in st.session_state:

            st.session_state["volcano_gene_text"] = ""

        if use_top_genes:

            top_n = st.number_input(
                "Top Up/Down Genes Per Direction",
                min_value=1,
                max_value=20,
                value=5,
                key="volcano_top_n"
            )

            #
            # Significant genes only
            #
            sig_df = de_df[
                de_df[significance_column]
                <= significance_cutoff
            ].copy()

            up_df_N = sum(sig_df["log2FC"]>= log2fc_cutoff)
            down_df_N = sum(sig_df["log2FC"]<= log2fc_cutoff)


            up_n = min(top_n,up_df_N)
            down_n = min(top_n,down_df_N)

            highlight_genes = sig_df['Gene'].tail(min(top_n,up_df_N)).to_list()+ sig_df['Gene'].head(min(top_n,down_df_N))

            #
            # Update only when changed
            #
            new_gene_text = ",".join(
                highlight_genes
            )

            if (
                st.session_state[
                    "volcano_gene_text"
                ]
                != new_gene_text
            ):

                st.session_state[
                    "highlight_genes"
                ] = highlight_genes

                st.session_state[
                    "volcano_gene_text"
                ] = new_gene_text

                st.session_state[
                    "boxplot_gene_text"
                ] = new_gene_text

                st.session_state[
                    "heatmap_gene_text"
                ] = new_gene_text

        #
        # Editable Gene List
        #
        gene_text = st.text_area(
            "Highlighted Genes",
            height=120,
            key="volcano_gene_text"
        )

        highlight_genes = [

            gene.strip()

            for gene in gene_text.split(",")

            if gene.strip()
        ]

        # --------------------------------------------------
        # Add Gene
        # --------------------------------------------------

        st.subheader(
            "Add Gene"
        )

        all_genes = sorted(
            de_df["Gene"]
            .astype(str)
            .unique()
            .tolist()
        )

        gene_to_add = st.selectbox(
            "Type Gene Name",
            options=all_genes,
            index=None,
            placeholder="Type to search...",
            key="volcano_gene_search"
        )

        if st.button(
            "Add Gene",
            key="volcano_add_gene"
        ):

            if (
                gene_to_add
                and gene_to_add
                not in highlight_genes
            ):

                highlight_genes = (
                    highlight_genes
                    + [gene_to_add]
                )

                new_text = ",".join(
                    highlight_genes
                )

                st.session_state[
                    "highlight_genes"
                ] = highlight_genes

                st.session_state[
                    "volcano_gene_text"
                ] = new_text

                st.session_state[
                    "boxplot_gene_text"
                ] = new_text

                st.session_state[
                    "heatmap_gene_text"
                ] = new_text
                st.write(highlight_genes)

                st.rerun()

        #
        # Keep synchronized
        #
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
                value=400,
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
                    step=0.1,
                    key="volcano_xmin"
                )

            with c6:
                x_max = st.number_input(
                    "X-axis Max",
                    value=5.0,
                    step=0.1,
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
                    step=0.1,
                    key="volcano_ymin"
                )

            with c8:
                y_max = st.number_input(
                    "Y-axis Max",
                    value=20.0,
                    step=0.1,
                    key="volcano_ymax"
                )

            y_range = [y_min, y_max]

        #
        # Create Figure
        #
        if st.button(
            "Generate Volcano Plot",
            key="generate_volcano_plot"
        ):
            st.write(highlight_genes)
            st.write(st.session_state.get("highlight_genes",[]))
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
    
    if (
        not st.session_state.get("boxplot_gene_text",[])
    ):

        st.session_state[
            "boxplot_gene_text"
        ] = ",".join(default_genes)
       

    gene_text = st.text_area(
        "Genes (comma separated)",
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
        meta_df.columns,
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
            step=50,
            key="boxplot_width"
        )

    with col2:

        plot_height = st.number_input(
            "Figure Height (px)",
            value=600,
            step=50,
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
                step=0.1,
                key="boxplot_ymin"
            )

        with col4:

            y_max = st.number_input(
                "Y-axis Maximum",
                value=20.0,
                step=0.1,
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
    if (
        not st.session_state.get("heatmap_genes",[])
    ):

        st.session_state[
            "heatmap_genes"
        ] = ",".join(default_genes)
        
    gene_text = st.text_area(
        "Genes (comma separated)",
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

    auto_color_scale = st.checkbox(
        "Automatic Color Scale",
        value=True,
        key="heatmap_auto_colorscale"
    )
    
    zmin = None
    zmid = None
    zmax = None

    if not auto_color_scale:

        c3, c4, c5 = st.columns(3)

        with c3:

            zmin = st.number_input(
                "Vmin",
                value=-2.0,
                step=0.1,
                key="heatmap_zmin"
            )

        with c4:

            zmid = st.number_input(
                "Center",
                value=0.0,
                step=0.1,
                key="heatmap_zmid"
            )

        with c5:

            zmax = st.number_input(
                "Vmax",
                value=2.0,
                step=0.1,
                key="heatmap_zmax"
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
        
    if st.button(
        "Generate Heatmap",
        key="generate_heatmap"
    ):
        if len(selected_genes) == 0:

            st.info(
                "Run DE analysis first or provide at least one gene."
            )

        else:
            fig = create_heatmap(
                expression_df=expr_df,
                metadata_df=meta_df,
                genes=selected_genes,
                annotation_column=annotation_column,
                apply_log2=apply_log2,
                zscore_by_gene=zscore_by_gene,
                cluster_samples=cluster_samples,
                cluster_genes=cluster_genes,
                width=heatmap_width,
                height=heatmap_height,
                colorscale=colorscale,
                zmin=zmin,
                zmid=zmid,
                zmax=zmax
            )
            import plotly
            if not isinstance(fig,plotly.graph_objs.Figure):
                st.write(fig)
                st.write(type(fig))
            else:

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
                
