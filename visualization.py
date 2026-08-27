# visualization.py

import numpy as np
import pandas as pd
#import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go

from scipy.cluster.hierarchy import (
    linkage,
    leaves_list
)
from scipy.stats import zscore

TAB10 = [
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # olive
    "#17becf",  # cyan
]


def apply_publication_style(
    fig
):
    if fig.layout.legend is not None:
        legend_gap_px = 40
        legend_x = (
            1
            + legend_gap_px / fig.layout.width
        )

        fig.update_layout(
            template="plotly_white",
            plot_bgcolor="white",
            paper_bgcolor="white",

            legend=dict(
                x=legend_x,
                y=1,
                xanchor="left",
                yanchor="top"
            )
        )
    else:
        fig.update_layout(
            template="plotly_white",
            plot_bgcolor="white",
            paper_bgcolor="white",
        )

    fig.update_xaxes(
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor="black"

    )

    fig.update_yaxes(
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        tickcolor="black"
      
    )
    return fig

def create_pca_plot(
    pca_df,
    color_column,
    explained_variance,
    width=500,
    height=500
):

    #
    # Build color map from categories
    #
    groups = (
        pca_df[color_column]
        .dropna()
        .unique()
    )

    color_map = {
        group: TAB10[i % len(TAB10)]
        for i, group in enumerate(groups)
    }

    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color=color_column,
        color_discrete_map=color_map,
        hover_data=["Sample"]
    )

    fig.update_traces(
        marker=dict(
            size=10
        )
    )

    fig.update_layout(
        template="plotly_white",
        width=width,
        height=height,
        xaxis_title=(
            f"PC1 ({explained_variance[0]*100:.1f}%)"
        ),
        yaxis_title=(
            f"PC2 ({explained_variance[1]*100:.1f}%)"
        ),
        legend_title=color_column
    )
    fig = apply_publication_style(fig)
    return fig


def create_volcano_plot(
    de_df,
    significance_column="FDR",
    significance_cutoff=0.05,
    log2fc_cutoff=1.0,
    highlight_genes=None,
    width=1200,
    height=800,
    x_range=None,
    y_range=None
):

    df = de_df.copy()

    df[significance_column] = (
        df[significance_column]
        .clip(lower=1e-300)
    )

    df["neglog10"] = (
        -np.log10(
            df[significance_column]
        )
    )

    df["Direction"] = "Not Significant"

    up_mask = (
        (df["log2FC"] >= log2fc_cutoff)
        &
        (
            df[significance_column]
            <= significance_cutoff
        )
    )

    down_mask = (
        (df["log2FC"] <= -log2fc_cutoff)
        &
        (
            df[significance_column]
            <= significance_cutoff
        )
    )

    df.loc[
        up_mask,
        "Direction"
    ] = "Up"

    df.loc[
        down_mask,
        "Direction"
    ] = "Down"

    color_map = {
        "Up": "red",
        "Down": "blue",
        "Not Significant": "lightgrey"
    }

    fig = px.scatter(
        df,
        x="log2FC",
        y="neglog10",
        color="Direction",
        color_discrete_map=color_map,
        hover_data={
            "Gene": True,
            "log2FC": ":.3f",
            "PValue": ":.3e",
            "FDR": ":.3e"
        }
    )

    fig.add_vline(
        x=log2fc_cutoff,
        line_dash="dash"
    )

    fig.add_vline(
        x=-log2fc_cutoff,
        line_dash="dash"
    )

    fig.add_hline(
        y=-np.log10(significance_cutoff),
        line_dash="dash"
    )

    if highlight_genes:

        highlight_df = df[
            df["Gene"].isin(
                highlight_genes
            )
        ]

        fig.add_trace(
            go.Scatter(
                x=highlight_df["log2FC"],
                y=highlight_df["neglog10"],
                mode="markers+text",
                text=highlight_df["Gene"],
                textposition="top center",
                marker=dict(
                    size=9,
                    color="purple",
                ),
                textfont=dict(
                    size=12,
                ),
                showlegend=False
            )
        )

    fig.update_layout(
        template="plotly_white",
        width=width,
        height=height,
        title="Volcano Plot",
        xaxis_title="log2 Fold Change",
        yaxis_title=f"-log10({significance_column})"
    )

    if x_range is not None:
        fig.update_xaxes(
            range=x_range
        )

    if y_range is not None:
        fig.update_yaxes(
            range=y_range
        )
    fig = apply_publication_style(fig)
    return fig
    
    
def create_gene_boxplot(
    expr_df,
    meta_df,
    group_column,
    width=1200,
    height=600,
    y_range=None,
    show_points=True,
    point_size=3,
    point_jitter=0.7
):
    """
    Create a gene-expression boxplot from pre-filtered data.

    Parameters
    ----------
    expr_df : pandas.DataFrame
        Pre-filtered expression matrix.

        Rows:
            Genes

        Columns:
            Samples

        The expression matrix must contain at least one gene
        and one sample.

    meta_df : pandas.DataFrame
        Pre-filtered metadata table.

        Rows:
            Samples

        The metadata index must match the expression-matrix
        columns.

    group_column : str
        Metadata column used to group and color samples.

    width : int
        Figure width in pixels.

    height : int
        Figure height in pixels.

    y_range : list or tuple or None
        Optional y-axis range:

            [y_min, y_max]

        If None, Plotly determines the range automatically.

    show_points : bool
        Whether to overlay individual sample points.

    point_size : int or float
        Size of the strip-plot points.

    point_jitter : float
        Horizontal jitter for the strip-plot points.

    Returns
    -------
    plotly.graph_objects.Figure
        Plotly boxplot figure.

    Notes
    -----
    Plot mode is inferred automatically:

        One gene:
            Single-gene boxplot

        More than one gene:
            Combined boxplot
    """

    # --------------------------------------------------
    # Convert expression data to long format
    # --------------------------------------------------

    plot_df = (
        expr_df
        .rename_axis("Gene")
        .reset_index()
        .melt(
            id_vars="Gene",
            var_name="Sample",
            value_name="Expression"
        )
    )

    plot_df = plot_df.merge(
        meta_df[[group_column]],
        left_on='Sample',
        right_index=True,
        how="inner",
        validate="many_to_one"
    )

    gene_order = (
        expr_df.index
        .tolist()
    )

    group_order = sorted(
        meta_df[group_column]
        .drop_duplicates()
        .tolist()
    )

    # --------------------------------------------------
    # Build TAB10 color mapping
    # --------------------------------------------------

    color_map = {
        group: TAB10[
            index % len(TAB10)
        ]
        for index, group
        in enumerate(group_order)
    }

    # --------------------------------------------------
    # Infer plot mode from number of genes
    # --------------------------------------------------

    single_gene_mode = (
        expr_df.shape[0] == 1
    )

    # ==================================================
    # SINGLE-GENE MODE
    # ==================================================

    if single_gene_mode:

        selected_gene = gene_order[0]

        fig = px.box(
            plot_df,
            x=group_column,
            y="Expression",
            color=group_column,
            color_discrete_map=color_map,
            category_orders={
                group_column: group_order
            },
            points=False,
            title=f"{selected_gene} Expression"
        )

        if show_points:

            strip_fig = px.strip(
                plot_df,
                x=group_column,
                y="Expression",
                color=group_column,
                category_orders={
                    group_column: group_order
                }
            )

            for trace in strip_fig.data:

                trace.marker.color = "black"
                trace.marker.size = point_size
                trace.marker.opacity = 0.7
                trace.jitter = point_jitter
                trace.showlegend = False
                trace.hovertemplate = (
                    "Expression=%{y}<extra></extra>"
                )

                fig.add_trace(
                    trace
                )

        x_axis_title = group_column

    # ==================================================
    # COMBINED MODE
    # ==================================================

    else:

        fig = px.box(
            plot_df,
            x="Gene",
            y="Expression",
            color=group_column,
            color_discrete_map=color_map,
            category_orders={
                "Gene": gene_order,
                group_column: group_order
            },
            points=False,
            title="Gene Expression"
        )

        if show_points:

            # Color must be supplied here even though all
            # markers are later changed to black. This makes
            # Plotly apply the same horizontal grouping offsets
            # used by the boxplot traces.
            strip_fig = px.strip(
                plot_df,
                x="Gene",
                y="Expression",
                color=group_column,
                category_orders={
                    "Gene": gene_order,
                    group_column: group_order
                }
            )

            for trace in strip_fig.data:

                trace.marker.color = "black"
                trace.marker.size = point_size
                trace.marker.opacity = 0.7
                trace.jitter = point_jitter
                trace.showlegend = False
                trace.hovertemplate = (
                    "Expression=%{y}<extra></extra>"
                )

                fig.add_trace(
                    trace
                )

        x_axis_title = "Gene"

    # --------------------------------------------------
    # Figure layout
    # --------------------------------------------------

    fig.update_layout(
        width=int(width),
        height=int(height),
        boxmode="group",
        xaxis_title=x_axis_title,
        yaxis_title="Expression",
        legend_title_text=group_column
    )

    if y_range is not None:

        if (
            len(y_range) != 2
            or y_range[0] >= y_range[1]
        ):

            raise ValueError(
                "y_range must contain a minimum and maximum, "
                "with minimum smaller than maximum."
            )

        fig.update_yaxes(
            range=[
                float(y_range[0]),
                float(y_range[1])
            ]
        )

    # Apply this last so that axes, ticks, backgrounds,
    # frames, and legend settings are consistent.
    fig = apply_publication_style(
        fig
    )

    return fig
    
def create_heatmap(
    expression_df,
    metadata_df,
    annotation_column,
    zscore_by_gene=True,
    cluster_samples=True,
    cluster_genes=True,
    width=1200,
    height=800,
    colorscale="RdBu_r",
    zmin=None,
    zmid=None,
    zmax=None,
):
    """
    Create an expression heatmap from pre-filtered data.

    Parameters
    ----------
    expression_df : pandas.DataFrame
        Filtered expression matrix.

        Rows:
            Genes

        Columns:
            Samples

    metadata_df : pandas.DataFrame
        Filtered metadata table.

        Rows:
            Samples

        The metadata index must match the expression
        matrix columns.

    annotation_column : str
        Metadata column displayed as a sample annotation
        bar above the expression heatmap.

    zscore_by_gene : bool
        If True, standardize each gene across samples.

    cluster_samples : bool
        If True, hierarchically cluster samples.

    cluster_genes : bool
        If True, hierarchically cluster genes.

    width : int
        Figure width in pixels.

    height : int
        Figure height in pixels.

    colorscale : str or list
        Plotly colorscale for expression values.

    zmin, zmid, zmax : float or None
        Optional expression color-scale limits and center.

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive expression heatmap with sample annotation.
    """

    # ==================================================
    # NUMERIC EXPRESSION MATRIX
    # ==================================================

    expr = expression_df.apply(
        pd.to_numeric,
        errors="coerce",
    )


    # ==================================================
    # GENE-WISE Z-SCORE
    # ==================================================

    if zscore_by_gene:

        row_means = expr.mean(
            axis=1,
        )

        row_stds = expr.std(
            axis=1,
            ddof=1,
        )

        row_stds = row_stds.replace(
            0,
            np.nan,
        )

        expr = (
            expr
            .sub(
                row_means,
                axis=0,
            )
            .div(
                row_stds,
                axis=0,
            )
            .fillna(0.0)
        )

    # ==================================================
    # CLUSTER GENES
    # ==================================================

    if (
        cluster_genes
        and expr.shape[0] > 1
    ):
        gene_linkage = linkage(
            expr.to_numpy(
                dtype=float,
            ),
            method="average",
            metric="euclidean",
        )

        gene_order = leaves_list(
            gene_linkage,
        )

        expr = expr.iloc[
            gene_order,
            :,
        ]

    # ==================================================
    # CLUSTER SAMPLES
    # ==================================================

    if (
        cluster_samples
        and expr.shape[1] > 1
    ):
        sample_linkage = linkage(
            expr.T.to_numpy(
                dtype=float,
            ),
            method="average",
            metric="euclidean",
        )

        sample_order = leaves_list(
            sample_linkage,
        )

        expr = expr.iloc[
            :,
            sample_order,
        ]

    # Reorder metadata to match the final sample order.
    metadata_df = metadata_df.loc[
        expr.columns
    ]

    # ==================================================
    # PREPARE SAMPLE ANNOTATION
    # ==================================================

    annotation_values = (
        metadata_df[
            annotation_column
        ]
    )

    annotation_is_numeric = (
        pd.api.types.is_numeric_dtype(
            annotation_values
        )
        and annotation_values.nunique(
            dropna=True
        ) > 10
    )

    if annotation_is_numeric:

        annotation_numeric = pd.to_numeric(
            annotation_values,
            errors="coerce",
        )

        annotation_fill_value = (
            annotation_numeric.median()
        )

        annotation_numeric = (
            annotation_numeric.fillna(
                annotation_fill_value,
            )
        )

        annotation_z = [
            annotation_numeric
            .to_numpy(
                dtype=float,
            )
            .tolist()
        ]

        annotation_colorscale = "Viridis"

        annotation_tickvals = None
        annotation_ticktext = None

        annotation_zmin = float(
            annotation_numeric.min()
        )

        annotation_zmax = float(
            annotation_numeric.max()
        )

        annotation_hover = np.array(
            [
                [
                    (
                        f"Sample: {sample}"
                        f"<br>{annotation_column}: "
                        f"{annotation_values.loc[sample]}"
                    )
                    for sample in expr.columns
                ]
            ]
        )

    else:

        annotation_text = (
            annotation_values
            .fillna("Missing")
            .astype(str)
        )

        annotation_categories = (
            annotation_text
            .drop_duplicates()
            .tolist()
        )

        annotation_code_map = {
            category: index
            for index, category
            in enumerate(
                annotation_categories
            )
        }

        annotation_codes = annotation_text.map(
            annotation_code_map
        )

        annotation_z = [
            annotation_codes
            .to_numpy(
                dtype=float,
            )
            .tolist()
        ]

        category_count = len(
            annotation_categories
        )

        if category_count == 1:
            annotation_zmin = -0.5
            annotation_zmax = 0.5
        else:
            annotation_zmin = -0.5
            annotation_zmax = (
                category_count - 0.5
            )

        # Construct a discrete Plotly colorscale
        # using the globally defined TAB10 palette.
        annotation_colorscale = []

        for index, category in enumerate(
            annotation_categories
        ):
            color = TAB10[
                index % len(TAB10)
            ]

            start = (
                index / category_count
            )

            end = (
                (index + 1)
                / category_count
            )

            annotation_colorscale.extend(
                [
                    [start, color],
                    [end, color],
                ]
            )

        annotation_tickvals = list(
            range(
                category_count
            )
        )

        annotation_ticktext = (
            annotation_categories
        )

        annotation_hover = np.array(
            [
                [
                    (
                        f"Sample: {sample}"
                        f"<br>{annotation_column}: "
                        f"{annotation_text.loc[sample]}"
                    )
                    for sample in expr.columns
                ]
            ]
        )

    # ==================================================
    # CREATE FIGURE
    # ==================================================

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.025,
        row_heights=[
            0.08,
            0.92,
        ],
    )

    # ==================================================
    # SAMPLE ANNOTATION BAR
    # ==================================================

    annotation_heatmap = go.Heatmap(
        z=annotation_z,
        x=expr.columns.tolist(),
        y=[annotation_column],
        colorscale=annotation_colorscale,
        zmin=annotation_zmin,
        zmax=annotation_zmax,
        customdata=annotation_hover,
        hovertemplate=(
            "%{customdata}"
            "<extra></extra>"
        ),
        showscale=True,
        colorbar={
            "title": {
                "text": annotation_column,
                "side": "right",
            },
            "x": 1.13,
            "y": 0.94,
            "len": 0.16,
            "thickness": 16,
            "outlinecolor": "black",
            "outlinewidth": 1,
            "tickvals": annotation_tickvals,
            "ticktext": annotation_ticktext,
        },
    )

    fig.add_trace(
        annotation_heatmap,
        row=1,
        col=1,
    )

    # ==================================================
    # EXPRESSION HEATMAP
    # ==================================================

    expression_colorbar_title = (
        "Gene Z-score"
        if zscore_by_gene
        else "Expression"
    )

    expression_heatmap = go.Heatmap(
        z=expr.to_numpy(
            dtype=float,
        ),
        x=expr.columns.tolist(),
        y=expr.index.tolist(),
        colorscale=colorscale,
        zmin=zmin,
        zmid=zmid,
        zmax=zmax,
        customdata=np.broadcast_to(
            expr.columns.to_numpy(),
            expr.shape,
        ),
        hovertemplate=(
            "Gene: %{y}"
            "<br>Sample: %{customdata}"
            "<br>Value: %{z:.3f}"
            "<extra></extra>"
        ),
        showscale=True,
        colorbar={
            "title": {
                "text": expression_colorbar_title,
                "side": "right",
            },
            "x": 1.02,
            "y": 0.44,
            "len": 0.72,
            "thickness": 16,
            "outlinecolor": "black",
            "outlinewidth": 1,
        },
    )

    fig.add_trace(
        expression_heatmap,
        row=2,
        col=1,
    )

    # ==================================================
    # FIGURE LAYOUT
    # ==================================================

    fig.update_layout(
        title={
            "text": "Gene Expression Heatmap",
            "x": 0.5,
            "xanchor": "center",
        },
        width=int(width),
        height=int(height),
        template="plotly_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin={
            "l": 120,
            "r": 240,
            "t": 80,
            "b": 140,
        },
    )

    # Annotation-axis formatting.
    fig.update_yaxes(
        title_text="",
        showgrid=False,
        showline=False,
        ticks="",
        row=1,
        col=1,
    )

    # Hide annotation-panel x labels because the same
    # sample labels are displayed under the main heatmap.
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        row=1,
        col=1,
    )

    # Expression-axis formatting.
    fig.update_xaxes(
        title_text="Samples",
        showgrid=False,
        ticks="outside",
        tickangle=45,
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        showline=True,
        linecolor="black",
        linewidth=1,
        mirror=True,
        row=2,
        col=1,
    )

    fig.update_yaxes(
        title_text="Genes",
        showgrid=False,
        ticks="outside",
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        showline=True,
        linecolor="black",
        linewidth=1,
        mirror=True,
        autorange="reversed",
        row=2,
        col=1,
    )

    return fig
    
def create_correlation_scatter(
    plot_df,
    x_column,
    y_column,
    method="pearson",
    statistics_df=None,
    group_column=None,
    width=700,
    height=600
):
    """
    Create a correlation scatter plot.

    If group_column is provided:
        - points are colored using TAB10
        - one fitted line is created per group
        - all selected groups appear together

    If group_column is None:
        - one overall scatter plot and fitted line
          are displayed
    """

    method = method.lower()

    trendline = (
        "ols"
        if method == "pearson"
        else "lowess"
    )

    plot_arguments = {
        "data_frame": plot_df,
        "x": x_column,
        "y": y_column,
        "trendline": trendline,
        "hover_data": ["Sample"]
    }

    # ----------------------------------
    # Group colors
    # ----------------------------------

    if (
        group_column is not None
        and group_column in plot_df.columns
    ):

        groups = (
            plot_df[group_column]
            .dropna()
            .unique()
            .tolist()
        )

        # Preserve the order represented in the data.
        color_map = {
            group: TAB10[
                index % len(TAB10)
            ]
            for index, group
            in enumerate(groups)
        }

        plot_arguments["color"] = (
            group_column
        )

        plot_arguments[
            "color_discrete_map"
        ] = color_map

        plot_arguments[
            "category_orders"
        ] = {
            group_column: groups
        }

    fig = px.scatter(
        **plot_arguments
    )

    # ----------------------------------
    # Title and statistics
    # ----------------------------------

    title = (
        f"{x_column} vs {y_column}"
    )

    if (
        statistics_df is not None
        and not statistics_df.empty
    ):

        statistic_name = (
            "r"
            if method == "pearson"
            else "rho"
        )

        statistic_lines = []

        for _, row in (
            statistics_df.iterrows()
        ):

            group_name = row.get(
                "Group",
                "All Samples"
            )

            coefficient = row.get(
                "Coefficient"
            )

            pvalue = row.get(
                "PValue"
            )

            fdr = row.get(
                "FDR"
            )

            n_samples = row.get(
                "N_Samples"
            )

            pieces = [
                f"{group_name}: "
                f"{statistic_name}="
                f"{coefficient:.3f}"
            ]

            if pd.notna(pvalue):

                pieces.append(
                    f"p={pvalue:.2e}"
                )

            if pd.notna(fdr):

                pieces.append(
                    f"FDR={fdr:.2e}"
                )

            if pd.notna(n_samples):

                pieces.append(
                    f"n={int(n_samples)}"
                )

            statistic_lines.append(
                " | ".join(pieces)
            )

        title += (
            "<br>"
            + "<br>".join(
                statistic_lines
            )
        )

    # ----------------------------------
    # Marker settings
    # ----------------------------------

    fig.update_traces(
        marker={
            "size": 8,
            "opacity": 0.8
        },
        selector={
            "mode": "markers"
        }
    )

    fig.update_layout(
        title={
            "text": title,
            "x": 0.5,
            "xanchor": "center"
        },
        width=width,
        height=height,
        legend_title_text=(
            group_column
            if group_column is not None
            else None
        ),
        margin={
            "l": 80,
            "r": 120,
            "t": (
                100
                if statistics_df is None
                else 130
                + 20 * len(
                    statistics_df
                )
            ),
            "b": 80
        }
    )

    fig = apply_publication_style(
        fig
    )

    return fig
    
# ==================================================
# Correlation Volcano Plot
# ==================================================

def create_correlation_volcano(
    corr_df,
    coefficient_cutoff=0.5,
    fdr_cutoff=0.05,
    width=1200,
    height=800
):

    df = corr_df.copy()

    df["minus_log10_fdr"] = (
        -np.log10(
            df["FDR"]
            .clip(lower=1e-300)
        )
    )

    df["Significant"] = (
        (
            df["Coefficient"]
            .abs()
            >= coefficient_cutoff
        )
        &
        (
            df["FDR"]
            <= fdr_cutoff
        )
    )

    fig = px.scatter(
        df,
        x="Coefficient",
        y="minus_log10_fdr",
        color="Significant",
        hover_data=[
            "Subject_A",
            "Subject_B",
            "Coefficient",
            "FDR"
        ]
    )

    fig.add_vline(
        x=coefficient_cutoff,
        line_dash="dash"
    )

    fig.add_vline(
        x=-coefficient_cutoff,
        line_dash="dash"
    )

    fig.add_hline(
        y=-np.log10(fdr_cutoff),
        line_dash="dash"
    )

    fig.update_layout(
        title="Correlation Volcano Plot",
        width=width,
        height=height,
        template="plotly_white"
    )
    fig = apply_publication_style(fig)
    return fig
    
# ==================================================
# Correlation Heatmap
# ==================================================

def create_correlation_heatmap(
    expr_df,
    genes,
    method="pearson",
    width=1000,
    height=900
):

    if len(genes) < 2:

        raise ValueError(
            "At least two genes are required."
        )

    gene_expr = expr_df.loc[
        genes
    ]

    corr_matrix = (
        gene_expr.T.corr(
            method=method
        )
    )

    fig = px.imshow(
        corr_matrix,
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        aspect="auto"
    )

    fig.update_layout(
        title=(
            f"{method.title()} "
            "Correlation Heatmap"
        ),
        width=width,
        height=height
    )
    fig = apply_publication_style(fig)
    return fig
    
# ==================================================
# Correlation Network Table
# ==================================================

def extract_network_edges(
    corr_df,
    coefficient_cutoff=0.7,
    fdr_cutoff=0.05
):

    pass    