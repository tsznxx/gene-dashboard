import streamlit as st
import pandas as pd
import numpy as np

from data_loader import (
    load_expression_file,
    load_metadata_file,
    validate_expression_matrix,
    validate_metadata,
    validate_sample_matching,
    summarize_expression,
    summarize_metadata,
)

from analysis import (
    apply_combat,
    run_pca,
    run_differential_expression,
    run_correlation_analysis,
    resolve_correlation_subject,
    get_correlation_plot_vectors
)

from visualization import (
    create_pca_plot,
    create_volcano_plot,
    create_gene_boxplot,
    create_heatmap,
    create_correlation_scatter,
    create_correlation_volcano,
    create_correlation_heatmap,
    #format_correlation_results
)
            
st.title("Gene Expression Dashboard")

st.set_page_config(
    page_title="Gene Expression Dashboard", page_icon="assets/logo.png", layout="wide"
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
    unsafe_allow_html=True,
)

# ==================================================
# SIDEBAR
# ==================================================

with st.sidebar:

    st.image("assets/logo.png", width="content")

    st.divider()
    st.header("Load Dataset")

    #
    # DATA ALREADY LOADED
    #
    if st.session_state.get("data_loaded", False):
        

        if st.button("Start Over", type="primary", key="start_over"):

            #
            # Clear entire session
            #
            for key in list(st.session_state.keys()):
                del st.session_state[key]

            st.rerun()

    #
    # NO DATA LOADED
    #
    else:
        # ----------------------
        # Example Data
        # ----------------------

        if st.button("Example Dataset", type="primary", key="load_example_dataset"):

            expr_df = pd.read_csv(
                "example_data/example_expression.tsv", sep="\t", index_col=0
            )
            meta_df = pd.read_csv(
                "example_data/example_metadata.tsv", sep="\t", index_col=0
            )

            st.session_state["expr_df"] = expr_df

            st.session_state["meta_df"] = meta_df

            st.session_state["data_loaded"] = True

            st.rerun()

        # ----------------------
        # Upload Data
        # ----------------------

        uploaded_expression = st.file_uploader(
            "Expression Matrix", type=["csv", "tsv", "txt"], key="expression_upload"
        )

        uploaded_metadata = st.file_uploader(
            "Metadata Table", type=["csv", "tsv", "txt"], key="metadata_upload"
        )

        if uploaded_expression is not None and uploaded_metadata is not None:

            expr_df = load_expression_file(uploaded_expression)

            meta_df = load_metadata_file(uploaded_metadata)

            st.session_state["expr_df"] = expr_df

            st.session_state["meta_df"] = meta_df

            st.session_state["data_loaded"] = True

            st.rerun()

    # ==================================================
    # Preprocess settings
    # ==================================================
    if st.session_state.get(
        "data_loaded",
        False
    ):

        st.divider()

        st.header(
            "Preprocessing"
        )

        # ----------------------------------
        # Initialize
        # ----------------------------------

        if "apply_log2" not in st.session_state:

            st.session_state[
                "apply_log2"
            ] = False

        if (
            "apply_batch_correction"
            not in st.session_state
        ):

            st.session_state[
                "apply_batch_correction"
            ] = False

        if (
            "current_preprocessing"
            not in st.session_state
        ):

            st.session_state[
                "current_preprocessing"
            ] = "raw"

        # ----------------------------------
        # User Requested Settings
        # ----------------------------------

        apply_log2 = st.checkbox(
            "Apply log2(x+1)",
            key="apply_log2"
        )

        apply_batch_correction = st.checkbox(
            "Apply Batch Correction (ComBat)",
            key="apply_batch_correction"
        )

        batch_column = None

        if apply_batch_correction:

            batch_candidates = [

                col

                for col in st.session_state[
                    "meta_df"
                ].columns

            ]

            #
            # meta_df already uses sample IDs
            # as index
            #
            if len(batch_candidates):

                default_index = 0

                for idx, col in enumerate(
                    batch_candidates
                ):

                    if (
                        col.lower()
                        == "batch"
                    ):

                        default_index = idx
                        break

                batch_column = st.selectbox(
                    "Batch Column",
                    batch_candidates,
                    index=default_index,
                    key="batch_column"
                )

        # ----------------------------------
        # Build requested key
        # ----------------------------------

        requested_steps = []

        if apply_log2:

            requested_steps.append(
                "log2"
            )

        if (
            apply_batch_correction
            and batch_column
        ):

            requested_steps.append(
                f"combat({batch_column})"
            )

        requested_key = (
            "->".join(
                requested_steps
            )
            if requested_steps
            else "raw"
        )

        # ----------------------------------
        # Current Status
        # ----------------------------------

        st.caption(
            "Current analysis matrix"
        )

        st.success(
            st.session_state[
                "current_preprocessing"
            ]
        )

        # ----------------------------------
        # Determine Whether Update Needed
        # ----------------------------------

        needs_update = (
            requested_key
            !=
            st.session_state[
                "current_preprocessing"
            ]
        )

        # ----------------------------------
        # Apply Button
        # ----------------------------------

        if needs_update:

            if st.button(
                "Apply Preprocessing",
                type="primary",
                key="apply_preprocessing"
            ):

                st.session_state[
                    "run_preprocessing"
                ] = True

                st.rerun()

        needs_update = (
            requested_key
            !=
            st.session_state[
                "current_preprocessing"
            ]
        )

    #
    # Example Files
    #
    st.divider()

    st.header("Example Files")

    with open("example_data/example_expression.tsv", "rb") as f:

        st.download_button(
            "Download Example Expression",
            data=f,
            file_name="example_expression.tsv",
            mime="text/tab-separated-values",
            key="download_example_expr",
        )

    with open("example_data/example_metadata.tsv", "rb") as f:

        st.download_button(
            "Download Example Metadata",
            data=f,
            file_name="example_metadata.tsv",
            mime="text/tab-separated-values",
            key="download_example_meta",
        )


#
# Main
#
# ==================================================
# LOAD DATA
# ==================================================

if not st.session_state.get(
    "data_loaded",
    False
):

    st.info(
        "Load the example dataset or upload your own dataset."
    )

    st.stop()


expr_df = st.session_state["expr_df"]
meta_df = st.session_state["meta_df"]

#
# Validate expression matrix
#
expr_validation = validate_expression_matrix(expr_df)


if expr_validation:

    st.error(expr_validation)

    st.stop()

#
# Validate metadata
#

meta_validation = validate_metadata(meta_df)


if meta_validation:

    st.error(meta_validation)

    st.stop()

#
# sample matching
#

matching_result = validate_sample_matching(expr_df, meta_df)

if not matching_result["matching"]:

    st.error("Sample mismatch detected.")

    if matching_result["missing_in_metadata"]:

        st.write("Missing in metadata:")

        st.write(matching_result["missing_in_metadata"])

    if matching_result["missing_in_expression"]:

        st.write("Missing in expression matrix:")

        st.write(matching_result["missing_in_expression"])

    st.stop()
meta_df = meta_df.reindex(index=expr_df.columns) # match the sample names.
st.success("Data loaded successfully.")


#
# Cache all genes
#
if "all_genes" not in st.session_state:

    st.session_state["all_genes"] = (
        expr_df.index.unique().astype(str)
        .tolist()
    )

# ==================================================
# PREPROCESSING SETUP
# ==================================================

#
# Store raw matrix once
#
if "expr_versions" not in st.session_state:

    st.session_state["expr_versions"] = {
        "raw": expr_df
    }

#
# First run
#
if "active_expr_df" not in st.session_state:

    st.session_state[
        "active_expr_df"
    ] = expr_df

if "current_preprocessing" not in st.session_state:

    st.session_state[
        "current_preprocessing"
    ] = "raw"

if "run_preprocessing" not in st.session_state:

    st.session_state[
        "run_preprocessing"
    ] = False

# ==================================================
# BUILD REQUESTED STATUS
# ==================================================

requested_steps = []

apply_log2 = st.session_state.get(
    "apply_log2",
    False
)

apply_batch_correction = st.session_state.get(
    "apply_batch_correction",
    False
)

batch_column = st.session_state.get(
    "batch_column"
)

if apply_log2:

    requested_steps.append(
        "log2"
    )

if (
    apply_batch_correction
    and batch_column
):

    requested_steps.append(
        f"combat({batch_column})"
    )

requested_key = (
    "->".join(
        requested_steps
    )
    if requested_steps
    else "raw"
)

# ==================================================
# RUN PREPROCESSING
# ==================================================

if st.session_state["run_preprocessing"]:

    #
    # Use cache if available
    #
    if (
        requested_key
        not in st.session_state[
            "expr_versions"
        ]
    ):

        with st.status(
            "Preprocessing dataset...",
            expanded=True
        ) as status:

            st.write(
                "Loading source matrix"
            )

            processed_df = st.session_state["expr_versions"]['raw']

            if apply_log2:

                st.write(
                    "Applying log2(x+1)"
                )
                if 'log2' in st.session_state["expr_versions"]:
                    processed_df = st.session_state["expr_versions"]['log2']
                else:
                    processed_df = np.log2(processed_df + 1)
                    st.session_state["expr_versions"]['log2'] = processed_df

            if apply_batch_correction:
                st.write(
                    f"Running ComBat: "
                    f"{batch_column}"
                )
                processed_df = apply_combat(
                    processed_df,
                    meta_df,
                    batch_column
                )

            st.write(
                "Saving result to cache"
            )

            status.update(
                label=(
                    f"Completed: "
                    f"{requested_key}"
                ),
                state="complete"
            )           
            #
            # Cache result
            #
            st.session_state[
                "expr_versions"
            ][requested_key] = processed_df
        
    #
    # Activate matrix
    #
    st.session_state[
        "active_expr_df"
    ] = st.session_state[
        "expr_versions"
    ][requested_key]

    #
    # Update current status
    #
    st.session_state[
        "current_preprocessing"
    ] = requested_key

    #
    # Reset flag
    #
    st.session_state[
        "run_preprocessing"
    ] = False
    st.rerun()

# ==================================================
# ACTIVE MATRIX FOR ANALYSIS
# ==================================================

expr_df = st.session_state[
    "active_expr_df"
]


#
# Tabs
#

tab_data, tab_pca, tab_correlation, tab_de_volcano, tab_box, tab_heatmap = st.tabs(
    ["Data", "PCA", "Correlation", "DEG Analysis", "Gene Boxplot", "Heatmap"]
)

#
# Data Tab
#

with tab_data:

    st.subheader("Dataset Summary")
    #st.write(st.session_state['expr_versions'].keys())

    expr_summary = summarize_expression(expr_df)

    meta_summary = summarize_metadata(meta_df)

    c1, c2, c3 = st.columns(3)

    with c1:
        st.metric("Genes", expr_summary["Genes"])

    with c2:
        st.metric("Samples", expr_summary["Samples"])

    with c3:
        st.metric("Groups", meta_summary.get("Groups", "NA"))

    st.divider()

    st.subheader("Expression Matrix Preview")

    st.dataframe(expr_df.head(20), width="content")

    st.subheader("Metadata Preview")

    st.dataframe(meta_df.head(20), width="content")

# ==================================================
# PCA TAB
# ==================================================

with tab_pca:

    st.subheader("Principal Component Analysis")

    color_column = st.selectbox(
        "Color samples by", meta_df.columns, key="pca_color_column"
    )

    st.divider()

    st.subheader("Figure Settings")

    col1, col2 = st.columns(2)

    with col1:

        pca_width = st.number_input(
            "Figure Width (px)",
            min_value=300,
            max_value=1200,
            value=500,
            step=50,
            key="pca_width",
        )

    with col2:

        pca_height = st.number_input(
            "Figure Height (px)",
            min_value=300,
            max_value=450,
            value=450,
            step=50,
            key="pca_height",
        )

    if st.button("Run PCA", key="run_pca_button"):

        try:

            pca_df, variance = run_pca(expr_df, apply_log2=apply_log2)

            pca_df = pca_df.merge(meta_df, on="Sample", how="left")

            fig = create_pca_plot(
                pca_df=pca_df,
                color_column=color_column,
                explained_variance=variance,
                width=pca_width,
                height=pca_height,
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
                        "scale": 1,
                    },
                },
            )

            st.subheader("Explained Variance")

            st.write(
                {
                    "PC1 (%)": round(variance[0] * 100, 2),
                    "PC2 (%)": round(variance[1] * 100, 2),
                }
            )

            st.subheader("PCA Coordinates")

            st.dataframe(pca_df, width="content")

        except Exception as e:

            st.error(str(e))


# ==================================================
# DIFFERENTIAL EXPRESSION AND VOLCANO TAB
# ==================================================

with tab_de_volcano:

    st.subheader(
        "Differential Expression and Volcano Plot"
    )

    current_preprocessing = (
        st.session_state.get(
            "current_preprocessing",
            "raw"
        )
    )

    st.caption(
        f"Current analysis matrix: "
        f"{current_preprocessing}"
    )

    # ==================================================
    # CHECK WHETHER STORED DE RESULTS ARE CURRENT
    # ==================================================

    de_results_current = (
        "de_results" in st.session_state
        and
        st.session_state.get(
            "de_preprocessing"
        )
        ==
        current_preprocessing
    )

    if (
        "de_results" in st.session_state
        and not de_results_current
    ):

        st.warning(
            "The stored differential expression results "
            "were generated from a different preprocessing "
            "state. Run DE analysis again for the current "
            "analysis matrix."
        )

    # ==================================================
    # DIFFERENTIAL EXPRESSION SETTINGS
    # ==================================================

    st.markdown(
        "### Differential Expression Settings"
    )

    # Exclude columns that cannot define two groups.
    candidate_group_columns = [
        column
        for column in meta_df.columns
        if (
            meta_df[column]
            .dropna()
            .nunique()
            >= 2
        )
    ]

    if len(candidate_group_columns) == 0:

        st.error(
            "No metadata column contains at least two "
            "groups. Differential expression cannot be run."
        )

    else:

        # Prefer a column named Group.
        group_default_index = 0

        for index, column in enumerate(
            candidate_group_columns
        ):

            if column.lower() == "group":

                group_default_index = index
                break

        group_column = st.selectbox(
            "Grouping Column",
            options=candidate_group_columns,
            index=group_default_index,
            key="de_group_column"
        )

        available_groups = (
            meta_df[group_column]
            .dropna()
            .unique()
            .tolist()
        )

        available_groups = sorted(
            available_groups,
            key=str
        )

        if len(available_groups) < 2:

            st.error(
                "The selected metadata column must contain "
                "at least two groups."
            )

        else:

            group_col1, group_col2 = st.columns(2)

            with group_col1:

                group1 = st.selectbox(
                    "Group 1",
                    options=available_groups,
                    index=0,
                    key="de_group1"
                )

            group2_options = [
                group
                for group in available_groups
                if group != group1
            ]

            with group_col2:

                group2 = st.selectbox(
                    "Group 2",
                    options=group2_options,
                    index=0,
                    key="de_group2"
                )

            group1_count = int(
                (
                    meta_df[group_column]
                    == group1
                ).sum()
            )

            group2_count = int(
                (
                    meta_df[group_column]
                    == group2
                ).sum()
            )

            st.caption(
                f"{group1}: {group1_count} samples | "
                f"{group2}: {group2_count} samples"
            )

            insufficient_samples = (
                group1_count < 2
                or group2_count < 2
            )

            if insufficient_samples:

                st.warning(
                    "Welch's t-test requires at least two "
                    "samples in each selected group."
                )

            # ==================================================
            # RUN DE ANALYSIS
            # ==================================================

            if st.button(
                "Run Differential Expression",
                type="primary",
                disabled=insufficient_samples,
                key="run_de_analysis"
            ):

                try:

                    with st.status(
                        "Running differential expression...",
                        expanded=True
                    ) as de_status:

                        de_status.write(
                            "Checking expression and metadata"
                        )

                        if expr_df is None:

                            raise ValueError(
                                "The active expression matrix "
                                "is unavailable."
                            )

                        if not isinstance(
                            expr_df,
                            pd.DataFrame
                        ):

                            raise TypeError(
                                "The active expression matrix "
                                "must be a pandas DataFrame."
                            )

                        if expr_df.empty:

                            raise ValueError(
                                "The active expression matrix "
                                "contains no genes."
                            )

                        if meta_df is None:

                            raise ValueError(
                                "The metadata table is unavailable."
                            )

                        if not isinstance(
                            meta_df,
                            pd.DataFrame
                        ):

                            raise TypeError(
                                "The metadata table must be "
                                "a pandas DataFrame."
                            )

                        if meta_df.empty:

                            raise ValueError(
                                "The metadata table is empty."
                            )

                        de_status.write(
                            f"Comparing {group1} versus {group2}"
                        )

                        de_results = (
                            run_differential_expression(
                                expression_df=expr_df,
                                metadata_df=meta_df,
                                group_column=group_column,
                                group1=group1,
                                group2=group2
                            )
                        )

                        if (
                            de_results is None
                            or de_results.empty
                        ):

                            raise ValueError(
                                "No valid differential expression "
                                "results were produced."
                            )

                        # Sort once and reuse the sorted table.
                        de_results = (
                            de_results
                            .sort_values(
                                "log2FC",
                                ascending=False
                            )
                            .reset_index(
                                drop=True
                            )
                        )

                        de_status.write(
                            "Saving differential expression results"
                        )

                        st.session_state[
                            "de_results"
                        ] = de_results

                        st.session_state[
                            "de_preprocessing"
                        ] = current_preprocessing

                        st.session_state[
                            "de_comparison"
                        ] = {
                            "group_column": group_column,
                            "group1": group1,
                            "group2": group2
                        }

                        # Reset Volcano top-gene signature so that
                        # defaults are regenerated for this analysis.
                        st.session_state.pop(
                            "volcano_top_signature",
                            None
                        )

                        st.session_state[
                            "highlight_genes"
                        ] = []

                        de_status.update(
                            label=(
                                "Differential expression completed"
                            ),
                            state="complete",
                            expanded=False
                        )

                    st.rerun()

                except Exception as error:

                    st.error(
                        "Differential expression analysis failed: "
                        f"{error}"
                    )

    # ==================================================
    # DE RESULTS
    # ==================================================

    de_results_current = (
        "de_results" in st.session_state
        and
        st.session_state.get(
            "de_preprocessing"
        )
        ==
        current_preprocessing
    )

    if de_results_current:

        de_results = st.session_state[
            "de_results"
        ]

        comparison = st.session_state.get(
            "de_comparison",
            {}
        )

        comparison_group1 = comparison.get(
            "group1",
            "Group 1"
        )

        comparison_group2 = comparison.get(
            "group2",
            "Group 2"
        )

        st.divider()

        st.markdown(
            "### Differential Expression Results"
        )

        st.write(
            f"**Comparison:** "
            f"{comparison_group1} vs "
            f"{comparison_group2}"
        )

        st.dataframe(
            de_results,
            width="stretch",
            hide_index=True
        )

        st.download_button(
            label="Download DE Results",
            data=de_results.to_csv(
                sep="\t",
                index=False
            ).encode(
                "utf-8"
            ),
            file_name=(
                f"DE_{comparison_group1}_vs_"
                f"{comparison_group2}.tsv"
            ),
            mime="text/tab-separated-values",
            key="download_de_results"
        )

        # ==================================================
        # VOLCANO SETTINGS
        # ==================================================

        st.divider()

        st.subheader("Volcano Plot")

        if not ("de_results" in st.session_state and st.session_state.get("de_preprocessing")==st.session_state.get("current_preprocessing")):

            st.info("Run Differential Expression first.")

        else:

            de_df = st.session_state["de_results"]
            if not "highlight_genes" in st.session_state:
                st.session_state["highlight_genes"] = []

            significance_column = st.radio(
                "Use significance metric",
                ["PValue", "FDR"],
                horizontal=True,
                key="volcano_sig_metric",
            )

            col1, col2 = st.columns(2)

            with col1:

                log2fc_cutoff = st.number_input(
                    "Absolute log2FC cutoff", value=1.0, step=0.1, key="volcano_fc_cutoff"
                )

            with col2:

                significance_cutoff = st.number_input(
                    f"{significance_column} cutoff",
                    value=0.05,
                    step=0.01,
                    key="volcano_sig_cutoff",
                )

            # --------------------------------------------------
            # Highlighted Genes
            # --------------------------------------------------

            st.subheader("Genes to Highlight")
            
            # --------------------------------------------------
            # Editable Gene List
            # --------------------------------------------------
            all_genes = st.session_state["all_genes"]
            gene_text = st.text_area(
                "Highlighted Genes",
                value=",".join(st.session_state["highlight_genes"]),
                height=120,
            )
            genes = [ gene.strip() for gene in gene_text.split(",") if gene.strip()]
            highlight_genes = [gene for gene in genes if gene in all_genes]
            st.session_state['not_found_genes'] = [gene for gene in genes if gene not in highlight_genes]
            st.session_state['highlight_genes'] = highlight_genes

            if "highlight_genes" not in st.session_state:
                st.session_state["highlight_genes"] = []

            use_top_genes = st.checkbox(
                "Use Top Significant Genes", value=True, key="volcano_use_top_genes"
            )

            if use_top_genes:

                top_n = st.number_input(
                    "Top Up/Down Genes Per Direction",
                    min_value=1,
                    max_value=100,
                    value=5,
                    key="volcano_top_n",
                )

                #
                # Significant genes only
                #
                sig_df = de_df[de_df[significance_column] <= significance_cutoff].copy()

                up_df_N = sum(sig_df["log2FC"] >= log2fc_cutoff)
                down_df_N = sum(sig_df["log2FC"] <= log2fc_cutoff)

                up_n = min(top_n, up_df_N)
                down_n = min(top_n, down_df_N)

                top_gene_list = (
                    sig_df["Gene"].tail(min(top_n, up_df_N)).to_list()
                    + sig_df["Gene"].head(min(top_n, down_df_N)).to_list()
                )

                top_gene_list = list(dict.fromkeys(top_gene_list))

                if st.button("Use Top Genes", key="refresh_top_genes"):

                    st.session_state["highlight_genes"] = top_gene_list

                    st.rerun()
                    
            # --------------------------------------------------
            # Add Gene
            # --------------------------------------------------

            colag1, colag2 = st.columns(2,vertical_alignment="bottom")
            with colag1:
                gene_to_add = st.selectbox(
                    "Type Gene Name",
                    options=[g for g in st.session_state.get("all_genes", []) if not g in st.session_state["highlight_genes"]] ,
                    index=None,
                    placeholder="Type to search...",
                    key="volcano_gene_search",
                )
            with colag2:
                if st.button("Add Gene", key="volcano_add_gene"):

                    if gene_to_add and gene_to_add not in highlight_genes:

                        st.session_state["highlight_genes"] = highlight_genes + [gene_to_add]

                        st.session_state["reload_volcano_text"] = True

                        st.rerun()
                    

            #
            # Keep synchronized
            #
            st.session_state["highlight_genes"] = highlight_genes
            not_found_genes = st.session_state.get('not_found_genes',[])
            if len(not_found_genes)>0:
                st.error(f'''Warning: [{",".join(not_found_genes)}] not found in genes!''')
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
                    key="volcano_width",
                )

            with col4:
                volcano_height = st.number_input(
                    "Figure Height (px)",
                    value=400,
                    min_value=200,
                    max_value=1500,
                    step=100,
                    key="volcano_height",
                )

            #
            # X-axis
            #
            auto_x = st.checkbox("Automatic X-axis", value=True, key="volcano_auto_x")

            x_range = None

            if not auto_x:

                c5, c6 = st.columns(2)

                with c5:
                    x_min = st.number_input(
                        "X-axis Min", value=-5.0, step=0.1, key="volcano_xmin"
                    )

                with c6:
                    x_max = st.number_input(
                        "X-axis Max", value=5.0, step=0.1, key="volcano_xmax"
                    )

                x_range = [x_min, x_max]

            #
            # Y-axis
            #
            auto_y = st.checkbox("Automatic Y-axis", value=True, key="volcano_auto_y")

            y_range = None

            if not auto_y:

                c7, c8 = st.columns(2)

                with c7:
                    y_min = st.number_input(
                        "Y-axis Min", value=0.0, step=0.1, key="volcano_ymin"
                    )

                with c8:
                    y_max = st.number_input(
                        "Y-axis Max", value=20.0, step=0.1, key="volcano_ymax"
                    )

                y_range = [y_min, y_max]

            #
            # Create Figure
            #
            if st.button("Generate Volcano Plot", key="generate_volcano_plot"):
                fig = create_volcano_plot(
                    de_df=de_df,
                    significance_column=significance_column,
                    significance_cutoff=significance_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    highlight_genes=highlight_genes,
                    width=volcano_width,
                    height=volcano_height,
                    x_range=x_range,
                    y_range=y_range,
                )

                #
                # Display Figure
                #
                st.plotly_chart(
                    fig,
                    width="content",
                    config={
                        "displaylogo": False,
                        "toImageButtonOptions": {
                            "format": "svg",
                            "filename": "volcano_plot",
                            "width": volcano_width,
                            "height": volcano_height,
                            "scale": 1,
                        },
                    },
                )

# ==================================================
# BOXPLOT TAB
# ==================================================

with tab_box:

    st.subheader("Gene Expression Boxplots")

    # ----------------------------------
    # Initialize
    # ----------------------------------

    if "boxplot_genes" not in st.session_state:

        st.session_state["boxplot_genes"] = st.session_state.get("highlight_genes", [])

    if "boxplot_gene_text" not in st.session_state:

        st.session_state["boxplot_gene_text"] = ",".join(
            st.session_state["boxplot_genes"]
        )

    if "reload_boxplot_text" not in st.session_state:

        st.session_state["reload_boxplot_text"] = False

    # ----------------------------------
    # Handle pending updates
    # ----------------------------------

    if st.session_state["reload_boxplot_text"]:

        st.session_state["boxplot_gene_text"] = ",".join(
            st.session_state["boxplot_genes"]
        )

        st.session_state["reload_boxplot_text"] = False
    # ----------------------------------
    # Text Area
    # ----------------------------------

    gene_text = st.text_area(
        "Genes (comma separated)", height=120, key="boxplot_gene_text"
    )

    selected_genes = [gene.strip() for gene in gene_text.split(",") if gene.strip() in st.session_state['all_genes']]

    # ----------------------------------
    # Buttons
    # ----------------------------------
    st.session_state["boxplot_genes"] = st.session_state.get("highlight_genes", [])

    if len(st.session_state.get("highlight_genes", []))>0:
        if st.button("Load Highlighted Genes", key="box_load_highlighted"):
            st.session_state["reload_boxplot_text"] = True
            st.rerun()

    # ----------------------------------
    # Add Gene
    # ----------------------------------

    colag1, colag2 = st.columns(2, vertical_alignment="bottom")
    with colag1:
        gene_to_add = st.selectbox(
            "Type Gene Name",
            options=st.session_state.get("all_genes", []),
            index=None,
            placeholder="Type to search...",
            key="boxplot_gene_search",
        )
    with colag2:
        if st.button("Add Gene", key="boxplot_add_gene"):

            if gene_to_add and gene_to_add not in selected_genes:

                st.session_state["boxplot_genes"] = selected_genes + [gene_to_add]

                st.session_state["reload_boxplot_text"] = True

                st.rerun()

    # ----------------------------------
    # Keep text edits
    # ----------------------------------

    st.session_state["boxplot_genes"] = selected_genes

    # ----------------------------------
    # Grouping
    # ----------------------------------

    group_column = st.selectbox(
        "Group By", meta_df.columns.tolist(), key="boxplot_group"
    )
    available_groups = (
        meta_df[group_column]
        .dropna()
        .unique()
        .tolist()
    )

    available_groups = sorted(
        available_groups,
        key=str
    )

    selected_groups = st.multiselect(
        "Groups to Include",
        options=available_groups,
        default=available_groups,
        key="corr_selected_groups"
    )

    if len(selected_groups) == 0:

        st.warning(
            "Select at least one group."
        )


    # ----------------------------------
    # Plot Mode
    # ----------------------------------

    plot_mode = st.radio("Plot Mode", ["Combined", "Single Gene"], horizontal=True)

    selected_gene = None

    if plot_mode == "Single Gene" and len(selected_genes) > 0:

        selected_gene = st.selectbox(
            "Select Gene", selected_genes, index=0, key="box_selected_gene"
        )

    # ----------------------------------
    # Figure Settings
    # ----------------------------------

    st.divider()

    st.subheader("Figure Settings")

    col3, col4 = st.columns(2)

    with col3:

        plot_width = st.number_input(
            "Figure Width (px)",
            value=1200,
            min_value=400,
            max_value=4000,
            step=100,
            key="boxplot_width",
        )

    with col4:

        plot_height = st.number_input(
            "Figure Height (px)",
            value=600,
            min_value=300,
            max_value=4000,
            step=100,
            key="boxplot_height",
        )

    auto_y = st.checkbox("Automatic Y-axis", value=True, key="boxplot_auto_y")

    y_range = None

    if not auto_y:

        col5, col6 = st.columns(2)

        with col5:

            y_min = st.number_input("Y-axis Minimum", value=0.0, key="boxplot_ymin")

        with col6:

            y_max = st.number_input("Y-axis Maximum", value=20.0, key="boxplot_ymax")

        y_range = [y_min, y_max]

    # ----------------------------------
    # Generate Plot
    # ----------------------------------

    if st.button("Generate Boxplot", key="generate_boxplot"):

        if len(selected_genes) == 0:

            st.warning("Please specify at least one gene.")

        else:
            genes = selected_genes if plot_mode == "Combined" else [selected_gene]
            selected_idx = meta_df[group_column].isin(selected_groups)
            filtered_expr_df = expr_df.loc[genes,selected_idx].copy()
            filtered_meta_df = meta_df.loc[selected_idx].copy()
            fig = create_gene_boxplot(
                expr_df=filtered_expr_df,
                meta_df=filtered_meta_df,
                group_column=group_column,
                width=plot_width,
                height=plot_height,
                y_range=y_range,
                show_points=True,
                point_size=3,
                point_jitter=0.7
            )

            filename = (
                "combined_boxplot"
                if plot_mode == "Combined"
                else f"{selected_gene}_boxplot"
            )

            st.plotly_chart(
                fig,
                width="content",
                config={
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "svg",
                        "filename": filename,
                        "width": plot_width,
                        "height": plot_height,
                        "scale": 1,
                    },
                },
            )

# ==================================================
# HEATMAP TAB
# ==================================================

with tab_heatmap:

    st.subheader("Expression Heatmap")

    # ----------------------------------
    # Initialize
    # ----------------------------------

    if "heatmap_genes" not in st.session_state:

        st.session_state["heatmap_genes"] = st.session_state.get("highlight_genes", [])

    if "heatmap_gene_text" not in st.session_state:

        st.session_state["heatmap_gene_text"] = ",".join(
            st.session_state["heatmap_genes"]
        )

    if "reload_heatmap_text" not in st.session_state:

        st.session_state["reload_heatmap_text"] = False

    # ----------------------------------
    # Handle pending updates
    # ----------------------------------

    if st.session_state["reload_heatmap_text"]:

        st.session_state["heatmap_gene_text"] = ",".join(
            st.session_state["heatmap_genes"]
        )

        st.session_state["reload_heatmap_text"] = False

    # ----------------------------------
    # Gene List
    # ----------------------------------

    gene_text = st.text_area(
        "Genes (comma separated)", height=120, key="heatmap_gene_text"
    )

    selected_genes = [gene.strip() for gene in gene_text.split(",") if gene.strip() in st.session_state['all_genes']]

    # ----------------------------------
    # Load highlighted genes
    # ----------------------------------
    if len(st.session_state.get("highlight_genes", []))>0:
        if st.button("Load Highlighted Genes", key="heatmap_load_highlighted"):

            st.session_state["heatmap_genes"] = st.session_state.get(
                "highlight_genes", []
            )

            st.session_state["reload_heatmap_text"] = True

            st.rerun()

    # ----------------------------------
    # Add Gene
    # ----------------------------------
    colag1, colag2 = st.columns(2, vertical_alignment="bottom")
    with colag1:
        gene_to_add = st.selectbox(
            "Type Gene Name",
            options=st.session_state.get("all_genes", []),
            index=None,
            placeholder="Type to search...",
            key="heatmap_gene_search",
        )
    with colag2:
        if st.button("Add Gene", key="heatmap_add_gene"):

            if gene_to_add and gene_to_add not in selected_genes:

                st.session_state["heatmap_genes"] = selected_genes + [gene_to_add]

                st.session_state["reload_heatmap_text"] = True

                st.rerun()

    # ----------------------------------
    # Save manual edits
    # ----------------------------------

    st.session_state["heatmap_genes"] = selected_genes

    # ----------------------------------
    # Heatmap Options
    # ----------------------------------

    annotation_column = st.selectbox(
        "Annotation Column", meta_df.columns.tolist(), index=1, key="heatmap_annotation"
    )

    col3, col4 = st.columns(2)

    with col3:

        zscore_by_gene = st.checkbox(
            "Z-score by Gene", value=True, key="heatmap_zscore"
        )

        cluster_genes = st.checkbox(
            "Cluster Genes", value=True, key="heatmap_cluster_genes"
        )

    with col4:

        cluster_samples = st.checkbox(
            "Cluster Samples", value=True, key="heatmap_cluster_samples"
        )

        colorscale = st.selectbox(
            "Color Map",
            ["RdBu_r", "Viridis", "Plasma", "Magma"],
            key="heatmap_colorscale",
        )

    # ----------------------------------
    # Color Scale Settings
    # ----------------------------------

    st.divider()

    st.subheader("Color Scale Settings")

    auto_color_scale = st.checkbox(
        "Automatic Color Scale", value=True, key="heatmap_auto_scale"
    )

    zmin = None
    zmid = None
    zmax = None

    if not auto_color_scale:

        c1, c2, c3 = st.columns(3)

        with c1:

            zmin = st.number_input("Vmin", value=-2.0, step=0.1, key="heatmap_zmin")

        with c2:

            zmid = st.number_input("Center", value=0.0, step=0.1, key="heatmap_zmid")

        with c3:

            zmax = st.number_input("Vmax", value=2.0, step=0.1, key="heatmap_zmax")

    # ----------------------------------
    # Figure Settings
    # ----------------------------------

    st.divider()

    st.subheader("Figure Settings")

    c4, c5 = st.columns(2)

    with c4:

        heatmap_width = st.number_input(
            "Figure Width (px)",
            value=1200,
            min_value=300,
            max_value=4000,
            step=100,
            key="heatmap_width",
        )

    with c5:

        heatmap_height = st.number_input(
            "Figure Height (px)",
            value=800,
            min_value=300,
            max_value=4000,
            step=100,
            key="heatmap_height",
        )

    # ----------------------------------
    # Generate Heatmap
    # ----------------------------------

    if st.button("Generate Heatmap", key="generate_heatmap"):

        selected_genes = st.session_state["heatmap_genes"]

        if len(selected_genes) == 0:

            st.warning("Please specify at least one gene.")

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
                zmax=zmax,
            )

            st.plotly_chart(
                fig,
                width="content",
                config={
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "svg",
                        "filename": "expression_heatmap",
                        "width": heatmap_width,
                        "height": heatmap_height,
                        "scale": 1,
                    },
                },
            )


# ==================================================
# CORRELATION TAB
# ==================================================

with tab_correlation:

    st.subheader(
        "Correlation Analysis"
    )

    st.caption(
        "The correlation analysis uses the current "
        f"analysis matrix: "
        f"{st.session_state.get('current_preprocessing', 'raw')}"
    )

    all_genes = sorted(
        st.session_state.get(
            "all_genes",
            expr_df.index.astype(str).tolist()
        )
    )

    subject_type_options = [
        "Single Gene",
        "Gene List",
        "Gene Signature",
        "All Genes"
    ]

    # ==================================================
    # SUBJECT TYPE SELECTION
    # ==================================================

    subject_type_col1, subject_type_col2 = (
        st.columns(2)
    )

    with subject_type_col1:

        subject_a_type = st.selectbox(
            "Subject A Type",
            options=subject_type_options,
            key="corr_subject_a_type"
        )

    with subject_type_col2:

        subject_b_type = st.selectbox(
            "Subject B Type",
            options=subject_type_options,
            key="corr_subject_b_type"
        )

    if (
        subject_a_type == "All Genes"
        and subject_b_type == "All Genes"
    ):

        st.warning(
            "All Genes versus All Genes is disabled "
            "because it can generate an extremely large "
            "number of pairwise comparisons. Select a "
            "Single Gene, Gene List, or Gene Signature "
            "for at least one subject."
        )

    st.divider()

    # ==================================================
    # SUBJECT A SETTINGS
    # ==================================================

    st.markdown(
        "### Subject A"
    )

    subject_a_gene = None
    subject_a_genes = None
    subject_a_signature_name = None
    subject_a_aggregation = "Mean"

    # --------------------------------------------------
    # Subject A: Single Gene
    # --------------------------------------------------

    if subject_a_type == "Single Gene":

        subject_a_gene = st.selectbox(
            "Subject A Gene",
            options=all_genes,
            index=0 if all_genes else None,
            key="corr_subject_a_gene"
        )

    # --------------------------------------------------
    # Subject A: Gene List
    # --------------------------------------------------

    elif subject_a_type == "Gene List":

        if "corr_a_gene_list" not in st.session_state:

            st.session_state[
                "corr_a_gene_list"
            ] = []

        if "corr_a_gene_text" not in st.session_state:

            st.session_state[
                "corr_a_gene_text"
            ] = ""

        if "corr_a_reload_text" not in st.session_state:

            st.session_state[
                "corr_a_reload_text"
            ] = False

        if st.session_state[
            "corr_a_reload_text"
        ]:

            st.session_state[
                "corr_a_gene_text"
            ] = ",".join(
                st.session_state[
                    "corr_a_gene_list"
                ]
            )

            st.session_state[
                "corr_a_reload_text"
            ] = False

        a_button_col1, a_button_col2 = (
            st.columns(2)
        )

        with a_button_col1:

            if st.button(
                "Load Highlighted Genes",
                key="corr_a_load_highlighted"
            ):

                st.session_state[
                    "corr_a_gene_list"
                ] = st.session_state.get(
                    "highlight_genes",
                    []
                ).copy()

                st.session_state[
                    "corr_a_reload_text"
                ] = True

                st.rerun()

        with a_button_col2:

            if st.button(
                "Clear Subject A Genes",
                key="corr_a_clear_genes"
            ):

                st.session_state[
                    "corr_a_gene_list"
                ] = []

                st.session_state[
                    "corr_a_reload_text"
                ] = True

                st.rerun()

        subject_a_gene_text = st.text_area(
            "Subject A Genes, comma separated",
            height=100,
            key="corr_a_gene_text"
        )

        subject_a_genes = list(
            dict.fromkeys(
                gene.strip()
                for gene
                in subject_a_gene_text.split(",")
                if gene.strip()
            )
        )

        subject_a_add_gene = st.selectbox(
            "Add a Gene to Subject A",
            options=all_genes,
            index=None,
            placeholder="Type to search for a gene...",
            key="corr_a_gene_search"
        )

        if st.button(
            "Add Gene to Subject A",
            key="corr_a_add_gene"
        ):

            if (
                subject_a_add_gene
                and subject_a_add_gene
                not in subject_a_genes
            ):

                st.session_state[
                    "corr_a_gene_list"
                ] = (
                    subject_a_genes
                    + [subject_a_add_gene]
                )

                st.session_state[
                    "corr_a_reload_text"
                ] = True

                st.rerun()

        st.session_state[
            "corr_a_gene_list"
        ] = subject_a_genes

    # --------------------------------------------------
    # Subject A: Gene Signature
    # --------------------------------------------------

    elif subject_a_type == "Gene Signature":

        subject_a_signature_name = st.text_input(
            "Subject A Signature Name",
            value="Signature A",
            key="corr_a_signature_name"
        )

        if (
            "corr_a_signature_genes"
            not in st.session_state
        ):

            st.session_state[
                "corr_a_signature_genes"
            ] = []

        if (
            "corr_a_signature_text"
            not in st.session_state
        ):

            st.session_state[
                "corr_a_signature_text"
            ] = ""

        if (
            "corr_a_signature_reload_text"
            not in st.session_state
        ):

            st.session_state[
                "corr_a_signature_reload_text"
            ] = False

        if st.session_state[
            "corr_a_signature_reload_text"
        ]:

            st.session_state[
                "corr_a_signature_text"
            ] = ",".join(
                st.session_state[
                    "corr_a_signature_genes"
                ]
            )

            st.session_state[
                "corr_a_signature_reload_text"
            ] = False

        a_sig_button_col1, a_sig_button_col2 = (
            st.columns(2)
        )

        with a_sig_button_col1:

            if st.button(
                "Load Highlighted Genes",
                key="corr_a_sig_load_highlighted"
            ):

                st.session_state[
                    "corr_a_signature_genes"
                ] = st.session_state.get(
                    "highlight_genes",
                    []
                ).copy()

                st.session_state[
                    "corr_a_signature_reload_text"
                ] = True

                st.rerun()

        with a_sig_button_col2:

            if st.button(
                "Clear Signature A Genes",
                key="corr_a_sig_clear"
            ):

                st.session_state[
                    "corr_a_signature_genes"
                ] = []

                st.session_state[
                    "corr_a_signature_reload_text"
                ] = True

                st.rerun()

        subject_a_signature_text = st.text_area(
            "Subject A Signature Genes, comma separated",
            height=100,
            key="corr_a_signature_text"
        )

        subject_a_genes = list(
            dict.fromkeys(
                gene.strip()
                for gene
                in subject_a_signature_text.split(",")
                if gene.strip()
            )
        )

        subject_a_signature_add_gene = st.selectbox(
            "Add a Gene to Signature A",
            options=all_genes,
            index=None,
            placeholder="Type to search for a gene...",
            key="corr_a_signature_gene_search"
        )

        if st.button(
            "Add Gene to Signature A",
            key="corr_a_signature_add_gene"
        ):

            if (
                subject_a_signature_add_gene
                and subject_a_signature_add_gene
                not in subject_a_genes
            ):

                st.session_state[
                    "corr_a_signature_genes"
                ] = (
                    subject_a_genes
                    + [
                        subject_a_signature_add_gene
                    ]
                )

                st.session_state[
                    "corr_a_signature_reload_text"
                ] = True

                st.rerun()

        st.session_state[
            "corr_a_signature_genes"
        ] = subject_a_genes

        subject_a_aggregation = st.selectbox(
            "Subject A Signature Scoring Method",
            options=[
                "Mean",
                "Median",
                "Sum",
                "Mean Z-score"
            ],
            key="corr_a_signature_aggregation"
        )

    elif subject_a_type == "All Genes":

        st.info(
            f"Subject A includes all "
            f"{len(all_genes):,} genes."
        )

    st.divider()

    # ==================================================
    # SUBJECT B SETTINGS
    # ==================================================

    st.markdown(
        "### Subject B"
    )

    subject_b_gene = None
    subject_b_genes = None
    subject_b_signature_name = None
    subject_b_aggregation = "Mean"

    # --------------------------------------------------
    # Subject B: Single Gene
    # --------------------------------------------------

    if subject_b_type == "Single Gene":

        default_b_index = (
            1 if len(all_genes) > 1 else 0
        )

        subject_b_gene = st.selectbox(
            "Subject B Gene",
            options=all_genes,
            index=(
                default_b_index
                if all_genes
                else None
            ),
            key="corr_subject_b_gene"
        )

    # --------------------------------------------------
    # Subject B: Gene List
    # --------------------------------------------------

    elif subject_b_type == "Gene List":

        if "corr_b_gene_list" not in st.session_state:

            st.session_state[
                "corr_b_gene_list"
            ] = []

        if "corr_b_gene_text" not in st.session_state:

            st.session_state[
                "corr_b_gene_text"
            ] = ""

        if "corr_b_reload_text" not in st.session_state:

            st.session_state[
                "corr_b_reload_text"
            ] = False

        if st.session_state[
            "corr_b_reload_text"
        ]:

            st.session_state[
                "corr_b_gene_text"
            ] = ",".join(
                st.session_state[
                    "corr_b_gene_list"
                ]
            )

            st.session_state[
                "corr_b_reload_text"
            ] = False

        b_button_col1, b_button_col2 = (
            st.columns(2)
        )

        with b_button_col1:

            if st.button(
                "Load Highlighted Genes",
                key="corr_b_load_highlighted"
            ):

                st.session_state[
                    "corr_b_gene_list"
                ] = st.session_state.get(
                    "highlight_genes",
                    []
                ).copy()

                st.session_state[
                    "corr_b_reload_text"
                ] = True

                st.rerun()

        with b_button_col2:

            if st.button(
                "Clear Subject B Genes",
                key="corr_b_clear_genes"
            ):

                st.session_state[
                    "corr_b_gene_list"
                ] = []

                st.session_state[
                    "corr_b_reload_text"
                ] = True

                st.rerun()

        subject_b_gene_text = st.text_area(
            "Subject B Genes, comma separated",
            height=100,
            key="corr_b_gene_text"
        )

        subject_b_genes = list(
            dict.fromkeys(
                gene.strip()
                for gene
                in subject_b_gene_text.split(",")
                if gene.strip()
            )
        )

        subject_b_add_gene = st.selectbox(
            "Add a Gene to Subject B",
            options=all_genes,
            index=None,
            placeholder="Type to search for a gene...",
            key="corr_b_gene_search"
        )

        if st.button(
            "Add Gene to Subject B",
            key="corr_b_add_gene"
        ):

            if (
                subject_b_add_gene
                and subject_b_add_gene
                not in subject_b_genes
            ):

                st.session_state[
                    "corr_b_gene_list"
                ] = (
                    subject_b_genes
                    + [subject_b_add_gene]
                )

                st.session_state[
                    "corr_b_reload_text"
                ] = True

                st.rerun()

        st.session_state[
            "corr_b_gene_list"
        ] = subject_b_genes

    # --------------------------------------------------
    # Subject B: Gene Signature
    # --------------------------------------------------

    elif subject_b_type == "Gene Signature":

        subject_b_signature_name = st.text_input(
            "Subject B Signature Name",
            value="Signature B",
            key="corr_b_signature_name"
        )

        if (
            "corr_b_signature_genes"
            not in st.session_state
        ):

            st.session_state[
                "corr_b_signature_genes"
            ] = []

        if (
            "corr_b_signature_text"
            not in st.session_state
        ):

            st.session_state[
                "corr_b_signature_text"
            ] = ""

        if (
            "corr_b_signature_reload_text"
            not in st.session_state
        ):

            st.session_state[
                "corr_b_signature_reload_text"
            ] = False

        if st.session_state[
            "corr_b_signature_reload_text"
        ]:

            st.session_state[
                "corr_b_signature_text"
            ] = ",".join(
                st.session_state[
                    "corr_b_signature_genes"
                ]
            )

            st.session_state[
                "corr_b_signature_reload_text"
            ] = False

        b_sig_button_col1, b_sig_button_col2 = (
            st.columns(2)
        )

        with b_sig_button_col1:

            if st.button(
                "Load Highlighted Genes",
                key="corr_b_sig_load_highlighted"
            ):

                st.session_state[
                    "corr_b_signature_genes"
                ] = st.session_state.get(
                    "highlight_genes",
                    []
                ).copy()

                st.session_state[
                    "corr_b_signature_reload_text"
                ] = True

                st.rerun()

        with b_sig_button_col2:

            if st.button(
                "Clear Signature B Genes",
                key="corr_b_sig_clear"
            ):

                st.session_state[
                    "corr_b_signature_genes"
                ] = []

                st.session_state[
                    "corr_b_signature_reload_text"
                ] = True

                st.rerun()

        subject_b_signature_text = st.text_area(
            "Subject B Signature Genes, comma separated",
            height=100,
            key="corr_b_signature_text"
        )

        subject_b_genes = list(
            dict.fromkeys(
                gene.strip()
                for gene
                in subject_b_signature_text.split(",")
                if gene.strip()
            )
        )

        subject_b_signature_add_gene = st.selectbox(
            "Add a Gene to Signature B",
            options=all_genes,
            index=None,
            placeholder="Type to search for a gene...",
            key="corr_b_signature_gene_search"
        )

        if st.button(
            "Add Gene to Signature B",
            key="corr_b_signature_add_gene"
        ):

            if (
                subject_b_signature_add_gene
                and subject_b_signature_add_gene
                not in subject_b_genes
            ):

                st.session_state[
                    "corr_b_signature_genes"
                ] = (
                    subject_b_genes
                    + [
                        subject_b_signature_add_gene
                    ]
                )

                st.session_state[
                    "corr_b_signature_reload_text"
                ] = True

                st.rerun()

        st.session_state[
            "corr_b_signature_genes"
        ] = subject_b_genes

        subject_b_aggregation = st.selectbox(
            "Subject B Signature Scoring Method",
            options=[
                "Mean",
                "Median",
                "Sum",
                "Mean Z-score"
            ],
            key="corr_b_signature_aggregation"
        )

    elif subject_b_type == "All Genes":

        st.info(
            f"Subject B includes all "
            f"{len(all_genes):,} genes."
        )

    st.divider()

    # ==================================================
    # ANALYSIS SETTINGS
    # ==================================================

    st.markdown(
        "### Analysis Settings"
    )

    run_by_group = st.checkbox(
        "Calculate correlations separately by group",
        value=False,
        key="corr_run_by_group"
    )

    group_column = None
    selected_groups = []

    if (run_by_group):

        eligible_group_columns = [
            column
            for column in meta_df.columns
            if (
                meta_df[column]
                .dropna()
                .nunique()
                >= 2
            )
        ]

        if len(eligible_group_columns) == 0:

            st.warning(
                "No metadata column contains at least two categories."
            )

        else:

            default_group_index = 0

            for index, column in enumerate(
                eligible_group_columns
            ):

                if column.lower() == "group":

                    default_group_index = index
                    break

            group_column = st.selectbox(
                "Group Category",
                options=eligible_group_columns,
                index=default_group_index,
                key="corr_group_column"
            )

            available_groups = (
                meta_df[group_column]
                .dropna()
                .unique()
                .tolist()
            )

            available_groups = sorted(
                available_groups,
                key=str
            )

            selected_groups = st.multiselect(
                "Groups to Include",
                options=available_groups,
                default=available_groups,
                key="corr_selected_groups"
            )

            if len(selected_groups) == 0:

                st.warning(
                    "Select at least one group."
                )

            group_counts = (
                meta_df.loc[
                    meta_df[group_column].isin(
                        selected_groups
                    ),
                    group_column
                ]
                .value_counts()
            )

            st.caption(
                " | ".join(
                    f"{group}: {count} samples"
                    for group, count
                    in group_counts.items()
                )
            )
            
    correlation_method = st.radio(
        "Correlation Method",
        options=[
            "pearson",
            "spearman"
        ],
        format_func=lambda value: value.title(),
        horizontal=True,
        key="corr_method"
    )

    min_samples = st.number_input(
        "Minimum Paired Samples",
        min_value=3,
        max_value=max(
            3,
            int(expr_df.shape[1])
        ),
        value=min(
            3,
            int(expr_df.shape[1])
        ),
        step=1,
        key="corr_min_samples"
    )

    skip_identical_pairs = st.checkbox(
        "Skip identical gene pairs",
        value=True,
        key="corr_skip_identical"
    )

    invalid_all_vs_all = (
        subject_a_type == "All Genes"
        and subject_b_type == "All Genes"
    )

    # ==================================================
    # RUN ANALYSIS
    # ==================================================

    if st.button(
        "Run Correlation Analysis",
        type="primary",
        disabled=invalid_all_vs_all,
        key="run_correlation_analysis"
    ):

        try:

            with st.status(
                "Running correlation analysis...",
                expanded=True
            ) as corr_status:

                corr_status.write(
                    "Resolving Subject A"
                )

                resolved_subject_a = (
                    resolve_correlation_subject(
                        expr_df=expr_df,
                        subject_type=
                        subject_a_type,
                        gene=subject_a_gene,
                        gene_list=
                        subject_a_genes,
                        signature_name=
                        subject_a_signature_name,
                        signature_aggregation=
                        subject_a_aggregation
                    )
                )

                corr_status.write(
                    "Resolving Subject B"
                )

                resolved_subject_b = (
                    resolve_correlation_subject(
                        expr_df=expr_df,
                        subject_type=
                        subject_b_type,
                        gene=subject_b_gene,
                        gene_list=
                        subject_b_genes,
                        signature_name=
                        subject_b_signature_name,
                        signature_aggregation=
                        subject_b_aggregation
                    )
                )

                if resolved_subject_a[
                    "missing_genes"
                ]:

                    st.warning(
                        "Subject A genes not found: "
                        + ", ".join(
                            resolved_subject_a[
                                "missing_genes"
                            ]
                        )
                    )

                if resolved_subject_b[
                    "missing_genes"
                ]:

                    st.warning(
                        "Subject B genes not found: "
                        + ", ".join(
                            resolved_subject_b[
                                "missing_genes"
                            ]
                        )
                    )

                corr_status.write(
                    "Computing pairwise correlations"
                )
                sample_groups = None

                if run_by_group:

                    if group_column is None:

                        raise ValueError(
                            "Select a group category before running "
                            "group-specific correlations."
                        )

                    if len(selected_groups) == 0:

                        raise ValueError(
                            "Select at least one group."
                        )

                    sample_groups = {}

                    for group_value in selected_groups:

                        samples_in_group = (
                            meta_df.index[
                                meta_df[group_column]
                                == group_value
                            ]
                            .astype(str)
                            .tolist()
                        )

                        sample_groups[
                            str(group_value)
                        ] = samples_in_group
        
                corr_results = run_correlation_analysis(
                    subject_a=resolved_subject_a,
                    subject_b=resolved_subject_b,
                    method=correlation_method,
                    min_samples=int(min_samples),
                    skip_identical_pairs=skip_identical_pairs,
                    allow_all_vs_all=False,
                    sample_groups=sample_groups
                )
                st.session_state["corr_group_settings"] = {
                    "run_by_group": run_by_group,
                    "group_column": group_column,
                    "selected_groups": (
                        list(selected_groups)
                        if selected_groups
                        else []
                    )
                }
                st.session_state[
                    "corr_results"
                ] = corr_results

                st.session_state[
                    "corr_subject_a"
                ] = resolved_subject_a

                st.session_state[
                    "corr_subject_b"
                ] = resolved_subject_b

                st.session_state[
                    "corr_method_used"
                ] = correlation_method

                st.session_state[
                    "corr_preprocessing"
                ] = st.session_state.get(
                    "current_preprocessing",
                    "raw"
                )

                corr_status.update(
                    label=(
                        "Correlation analysis completed"
                    ),
                    state="complete",
                    expanded=False
                )

        except Exception as error:

            st.error(
                f"Correlation analysis failed: "
                f"{error}"
            )

    # ==================================================
    # CHECK WHETHER RESULTS ARE CURRENT
    # ==================================================

    corr_results_current = (
        "corr_results" in st.session_state
        and
        st.session_state.get(
            "corr_preprocessing"
        )
        ==
        st.session_state.get(
            "current_preprocessing",
            "raw"
        )
    )

    if (
        "corr_results" in st.session_state
        and not corr_results_current
    ):

        st.warning(
            "The stored correlation results were "
            "generated from a different preprocessing "
            "state. Run correlation analysis again."
        )

    # ==================================================
    # RESULTS
    # ==================================================

    if corr_results_current:

        corr_results = st.session_state[
            "corr_results"
        ]

        st.divider()

        st.markdown(
            "### Correlation Results"
        )

        if corr_results.empty:

            st.warning(
                "No valid correlation results were produced."
            )

        else:

            result_col1, result_col2 = (
                st.columns(2)
            )

            with result_col1:

                maximum_rows = max(
                    1,
                    len(corr_results)
                )

                default_rows = min(
                    100,
                    maximum_rows
                )
                rows_to_show = st.number_input(
                    f"Rows to Display (out of {len(corr_results)})",
                    min_value=1,
                    max_value=maximum_rows,
                    value=default_rows,
                    step=1,
                    key="corr_rows_to_show"
                )

            with result_col2:

                absolute_cutoff = st.number_input(
                    "Minimum Absolute Correlation",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.0,
                    step=0.05,
                    key="corr_result_cutoff"
                )

            filtered_corr_results = (
                corr_results[
                    corr_results[
                        "Absolute_Correlation"
                    ]
                    >= absolute_cutoff
                ]
                .head(
                    int(rows_to_show)
                )
            )

            display_corr_results = (
                filtered_corr_results.copy()
            )

            display_corr_results[
                "Coefficient"
            ] = display_corr_results[
                "Coefficient"
            ].round(4)

            display_corr_results[
                "PValue"
            ] = display_corr_results[
                "PValue"
            ].map(
                lambda value:
                (
                    f"{value:.3e}"
                    if pd.notna(value)
                    else ""
                )
            )

            display_corr_results[
                "FDR"
            ] = display_corr_results[
                "FDR"
            ].map(
                lambda value:
                (
                    f"{value:.3e}"
                    if pd.notna(value)
                    else ""
                )
            )

            st.dataframe(
                display_corr_results[
                    [
                        "Group",
                        "Subject_A",
                        "Subject_B",
                        "Coefficient",
                        "PValue",
                        "FDR",
                        "N_Samples"
                    ]
                ],
                width="stretch",
                hide_index=True
            )

            st.download_button(
                label="Download Correlation Results",
                data=corr_results.to_csv(
                    sep="\t",
                    index=False
                ).encode(
                    "utf-8"
                ),
                file_name=(
                    "correlation_results.tsv"
                ),
                mime=(
                    "text/tab-separated-values"
                ),
                key="download_corr_results"
            )

            # ==================================================
            # CORRELATION SCATTER PLOT
            # ==================================================

            st.divider()

            st.markdown(
                "### Correlation Scatter Plot"
            )

            # --------------------------------------------------
            # Select Subject A
            # --------------------------------------------------

            subject_a_options = (
                corr_results["Subject_A"]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            if len(subject_a_options) == 0:

                st.warning(
                    "No Subject A values are available "
                    "for scatter plotting."
                )

            else:

                selected_scatter_a = st.selectbox(
                    "Scatter Plot Subject A",
                    options=subject_a_options,
                    key="corr_scatter_subject_a"
                )

                # --------------------------------------------------
                # Select Subject B
                # --------------------------------------------------

                subject_b_options = (
                    corr_results.loc[
                        corr_results[
                            "Subject_A"
                        ].astype(str)
                        == str(
                            selected_scatter_a
                        ),
                        "Subject_B"
                    ]
                    .dropna()
                    .astype(str)
                    .unique()
                    .tolist()
                )

                if len(subject_b_options) == 0:

                    st.warning(
                        "No Subject B values are available "
                        "for the selected Subject A."
                    )

                else:

                    selected_scatter_b = st.selectbox(
                        "Scatter Plot Subject B",
                        options=subject_b_options,
                        key="corr_scatter_subject_b"
                    )

                    # ==================================================
                    # FIGURE SETTINGS
                    # ==================================================

                    st.markdown(
                        "#### Figure Settings"
                    )

                    figure_col1, figure_col2 = (
                        st.columns(2)
                    )

                    with figure_col1:

                        corr_plot_width = st.number_input(
                            "Scatter Plot Width (px)",
                            min_value=300,
                            max_value=2000,
                            value=700,
                            step=50,
                            key="corr_plot_width"
                        )

                    with figure_col2:

                        corr_plot_height = st.number_input(
                            "Scatter Plot Height (px)",
                            min_value=300,
                            max_value=2000,
                            value=600,
                            step=50,
                            key="corr_plot_height"
                        )

                    # ==================================================
                    # RETRIEVE SAVED GROUP SETTINGS
                    # ==================================================

                    group_settings = (
                        st.session_state.get(
                            "corr_group_settings",
                            {}
                        )
                    )

                    stored_run_by_group = (
                        group_settings.get(
                            "run_by_group",
                            False
                        )
                    )

                    stored_group_column = (
                        group_settings.get(
                            "group_column"
                        )
                    )

                    stored_selected_groups = (
                        group_settings.get(
                            "selected_groups",
                            []
                        )
                    )

                    # ==================================================
                    # RETRIEVE EXPRESSION VECTORS
                    # ==================================================

                    try:

                        scatter_data = (
                            get_correlation_plot_vectors(
                                subject_a=
                                st.session_state[
                                    "corr_subject_a"
                                ],
                                subject_b=
                                st.session_state[
                                    "corr_subject_b"
                                ],
                                subject_a_name=
                                selected_scatter_a,
                                subject_b_name=
                                selected_scatter_b,
                                sample_ids=None
                            )
                        )

                        # Normalize sample IDs before
                        # merging with metadata.
                        scatter_data["Sample"] = (
                            scatter_data["Sample"]
                            .astype(str)
                        )

                        # ==================================================
                        # ADD GROUP METADATA
                        # ==================================================

                        scatter_group_column = None

                        if (
                            stored_run_by_group
                            and stored_group_column is not None
                            and stored_group_column
                            in meta_df.columns
                        ):

                            scatter_group_column = (
                                stored_group_column
                            )

                            annotation_df = (
                                meta_df[
                                    [
                                        stored_group_column
                                    ]
                                ]
                                .copy()
                            )

                            annotation_df.index = (
                                annotation_df.index
                                .astype(str)
                            )

                            annotation_df.index.name = (
                                "Sample"
                            )

                            annotation_df = (
                                annotation_df
                                .reset_index()
                            )

                            scatter_data = (
                                scatter_data.merge(
                                    annotation_df,
                                    on="Sample",
                                    how="left"
                                )
                            )

                            # Keep only the groups that were selected
                            # when the correlation analysis was run.
                            if (
                                stored_selected_groups
                                is not None
                                and len(
                                    stored_selected_groups
                                ) > 0
                            ):

                                scatter_data = (
                                    scatter_data[
                                        scatter_data[
                                            stored_group_column
                                        ].isin(
                                            stored_selected_groups
                                        )
                                    ]
                                    .copy()
                                )

                        # ==================================================
                        # RETRIEVE STATISTICS FOR THE SELECTED PAIR
                        # ==================================================

                        scatter_statistics = (
                            corr_results[
                                (
                                    corr_results[
                                        "Subject_A"
                                    ].astype(str)
                                    == str(
                                        selected_scatter_a
                                    )
                                )
                                &
                                (
                                    corr_results[
                                        "Subject_B"
                                    ].astype(str)
                                    == str(
                                        selected_scatter_b
                                    )
                                )
                            ]
                            .copy()
                        )

                        # Preserve the group order chosen by the
                        # user when group-specific analysis was run.
                        if (
                            stored_run_by_group
                            and stored_selected_groups
                            and "Group"
                            in scatter_statistics.columns
                        ):

                            group_order = {
                                str(group):
                                index
                                for index, group
                                in enumerate(
                                    stored_selected_groups
                                )
                            }

                            scatter_statistics[
                                "_group_order"
                            ] = (
                                scatter_statistics["Group"]
                                .astype(str)
                                .map(group_order)
                            )

                            scatter_statistics = (
                                scatter_statistics
                                .sort_values(
                                    "_group_order"
                                )
                                .drop(
                                    columns=[
                                        "_group_order"
                                    ]
                                )
                            )

                        # ==================================================
                        # VALIDATE PLOTTING DATA
                        # ==================================================

                        if scatter_data.empty:

                            st.warning(
                                "No samples are available for "
                                "the selected correlation pair."
                            )

                        elif (
                            selected_scatter_a
                            not in scatter_data.columns
                            or selected_scatter_b
                            not in scatter_data.columns
                        ):

                            st.warning(
                                "The selected correlation vectors "
                                "are unavailable in the plotting data."
                            )

                        else:

                            # Display sample counts for grouped plots.
                            if (
                                scatter_group_column is not None
                                and scatter_group_column
                                in scatter_data.columns
                            ):

                                scatter_group_counts = (
                                    scatter_data[
                                        scatter_group_column
                                    ]
                                    .value_counts()
                                )

                                st.caption(
                                    " | ".join(
                                        (
                                            f"{group}: "
                                            f"{count} samples"
                                        )
                                        for group, count
                                        in scatter_group_counts.items()
                                    )
                                )

                            # ==================================================
                            # CREATE SCATTER FIGURE
                            # ==================================================
                            if st.button("Generate Scatter Plot", key="generate_correlation_scatter"):
                                scatter_fig = (
                                    create_correlation_scatter(
                                        plot_df=scatter_data,
                                        x_column=
                                        selected_scatter_a,
                                        y_column=
                                        selected_scatter_b,
                                        method=
                                        st.session_state.get(
                                            "corr_method_used",
                                            "pearson"
                                        ),
                                        statistics_df=
                                        scatter_statistics,
                                        group_column=
                                        scatter_group_column,
                                        width=int(
                                            corr_plot_width
                                        ),
                                        height=int(
                                            corr_plot_height
                                        )
                                    )
                                )

                                # ==================================================
                                # EXPORT FILENAME
                                # ==================================================

                                scatter_filename = (
                                    f"{selected_scatter_a}_vs_"
                                    f"{selected_scatter_b}_"
                                    f"correlation"
                                )

                                if (
                                    stored_run_by_group
                                    and stored_group_column
                                ):

                                    scatter_filename += (
                                        f"_by_{stored_group_column}"
                                    )

                                # Replace filename characters that may
                                # be problematic on some systems.
                                scatter_filename = (
                                    scatter_filename
                                    .replace(
                                        " ",
                                        "_"
                                    )
                                    .replace(
                                        "/",
                                        "_"
                                    )
                                    .replace(
                                        "\\",
                                        "_"
                                    )
                                )

                                # ==================================================
                                # DISPLAY AND EXPORT SVG
                                # ==================================================

                                st.plotly_chart(
                                    scatter_fig,
                                    width="content",
                                    config={
                                        "displaylogo": False,
                                        "toImageButtonOptions": {
                                            "format": "svg",
                                            "filename":
                                            scatter_filename,
                                            "width":
                                            int(
                                                corr_plot_width
                                            ),
                                            "height":
                                            int(
                                                corr_plot_height
                                            ),
                                            "scale": 1
                                        }
                                    }
                                )

                    except Exception as error:

                        st.error(
                            "Unable to create the correlation "
                            f"scatter plot: {error}"
                        )