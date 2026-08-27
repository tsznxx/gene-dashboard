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
    expression_df,
    metadata_df,
    genes,
    group_column,
    apply_log2=False,
    plot_mode="Combined",
    selected_gene=None,
    width=1200,
    height=600,
    y_range=None
):

    expr = expression_df

    expr = expr.apply(
        pd.to_numeric,
        errors="coerce"
    )

    plot_records = []

    for gene in genes:

        if gene not in expr.index:
            continue

        for sample in expr.columns:

            plot_records.append(
                {
                    "Gene": gene,
                    "Sample": sample,
                    "Expression": expr.loc[
                        gene,
                        sample
                    ]
                }
            )

    plot_df = pd.DataFrame(
        plot_records
    )

    plot_df = plot_df.merge(
        metadata_df,
        left_index=True,
        right_index=True,
        how="left"
    )

    groups = (
        plot_df[group_column]
        .dropna()
        .unique()
    )

    color_map = {
        group: TAB10[i % len(TAB10)]
        for i, group in enumerate(groups)
    }

    if plot_mode == "Single Gene":

        if selected_gene is None:
            selected_gene = genes[0]

        plot_df = plot_df[
            plot_df["Gene"]
            ==
            selected_gene
        ]

        fig = px.box(
            plot_df,
            x=group_column,
            y="Expression",
            color=group_column,
            color_discrete_map=color_map,
            points=False,
            title=f"{selected_gene} Expression"
        )

        strip_fig = px.strip(
            plot_df,
            x=group_column,
            y="Expression"
        )

        for trace in strip_fig.data:

            trace.marker.color = "black"
            trace.marker.size = 3
            trace.jitter = 1.0
            trace.showlegend = False

            fig.add_trace(trace)

    else:
        fig = px.box(
            plot_df,
            x="Gene",
            y="Expression",
            color=group_column,
            color_discrete_map=color_map,
            points=False
        )

        strip_fig = px.strip(
            plot_df,
            x="Gene",
            y="Expression",
            color=group_column
        )

        for trace in strip_fig.data:

            trace.marker.color = "black"
            trace.marker.size = 3
            trace.jitter = 0.8
            trace.showlegend = False

            fig.add_trace(trace)


    fig.update_layout(
        template="plotly_white",
        width=width,
        height=height
    )

    if y_range is not None:

        fig.update_yaxes(
            range=y_range
        )
    fig = apply_publication_style(fig)
    return fig
    
def create_heatmap(
    expression_df,
    metadata_df,
    genes,
    annotation_column,
    apply_log2=False,
    zscore_by_gene=True,
    cluster_samples=True,
    cluster_genes=True,
    width=1200,
    height=800,
    colorscale="RdBu_r",
    zmin=None,
    zmid=None,
    zmax=None
):


    #
    # Expression matrix
    #

    expr = expression_df

    expr = expr.apply(
        pd.to_numeric,
        errors="coerce"
    )

    #
    # Keep selected genes
    #
    valid_genes = [
        g for g in genes
        if g in expr.index
    ]

    expr = expr.loc[
        valid_genes
    ]

    #
    # Z-score by gene
    #
    if zscore_by_gene:

        expr = expr.sub(
            expr.mean(axis=1),
            axis=0
        ).div(
            expr.std(axis=1),
            axis=0
        )

        expr = expr.fillna(0)
    #
    # Cluster genes
    #
    if (
        cluster_genes
        and expr.shape[0] > 1
    ):

        gene_linkage = linkage(
            expr,
            method="average"
        )

        gene_order = leaves_list(
            gene_linkage
        )

        expr = expr.iloc[
            gene_order
        ]

    #
    # Cluster samples
    #
    if (
        cluster_samples
        and expr.shape[1] > 1
    ):

        sample_linkage = linkage(
            expr.T,
            method="average"
        )

        sample_order = leaves_list(
            sample_linkage
        )

        expr = expr.iloc[
            :,
            sample_order
        ]

    #
    # Sample labels
    #
    sample_labels = []

    annotation_map = (
        metadata_df
        [annotation_column]
        .to_dict()
    )

    for sample in expr.columns:

        value = annotation_map.get(
            sample,
            ""
        )

        sample_labels.append(
            f"{sample}<br>{value}"
        )

    fig = px.imshow(
        expr,
        color_continuous_scale=colorscale,
        zmin=zmin,
        zmax=zmax,
        aspect="auto"
    )
    if zmid is not None:

        fig.update_coloraxes(
            cmid=zmid
        )

    fig.update_xaxes(
        tickvals=list(
            range(
                len(expr.columns)
            )
        ),
        ticktext=sample_labels
    )

    fig.update_layout(
        width=width,
        height=height,
        template="plotly_white",
        xaxis_title="Samples",
        yaxis_title="Genes"
    )

    #fig = apply_publication_style(fig)
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