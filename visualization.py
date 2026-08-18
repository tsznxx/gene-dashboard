# visualization.py

import numpy as np
import pandas as pd
#import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go


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


def create_pca_plot(
    pca_df,
    color_column,
    explained_variance
):
    """
    Generate interactive PCA plot.
    """

    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color=color_column,
        hover_data=["Sample"],
        title="Principal Component Analysis"
    )

    fig.update_traces(
        marker=dict(
            size=10
        )
    )

    fig.update_layout(
        template="plotly_white",
        height=700,
        xaxis_title=(
            f"PC1 "
            f"({explained_variance[0]*100:.1f}%)"
        ),
        yaxis_title=(
            f"PC2 "
            f"({explained_variance[1]*100:.1f}%)"
        )
    )

    return fig


def create_volcano_plot(
    de_df,
    significance_column="FDR",
    significance_cutoff=0.05,
    log2fc_cutoff=1.0,
    highlight_genes=None
):
    """
    Create volcano plot.

    Parameters
    ----------
    de_df : DataFrame

        Must contain:

        Gene
        log2FC
        PValue
        FDR

    significance_column : str

        "PValue" or "FDR"

    significance_cutoff : float

    log2fc_cutoff : float

    highlight_genes : list
    """

    df = de_df.copy()

    #
    # Avoid log(0)
    #
    df[significance_column] = (
        df[significance_column]
        .clip(lower=1e-300)
    )

    df["neglog10"] = (
        -np.log10(
            df[significance_column]
        )
    )

    #
    # Define Up/Down/NS
    #
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
            "FDR": ":.3e",
            "Direction": True
        }
    )

    #
    # FC thresholds
    #
    fig.add_vline(
        x=log2fc_cutoff,
        line_dash="dash",
        line_color="black"
    )

    fig.add_vline(
        x=-log2fc_cutoff,
        line_dash="dash",
        line_color="black"
    )

    #
    # Significance threshold
    #
    fig.add_hline(
        y=-np.log10(significance_cutoff),
        line_dash="dash",
        line_color="black"
    )

    #
    # Highlight genes
    #
    if (
        highlight_genes is not None
        and len(highlight_genes) > 0
    ):

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
                    color="purple",
                    size=11,
                    line=dict(
                        color="white",
                        width=1
                    )
                ),
                showlegend=False,
                hoverinfo="skip"
            )
        )

    fig.update_layout(
        template="plotly_white",
        height=750,
        title=(
            f"Volcano Plot "
            f"({significance_column} ≤ {significance_cutoff}, "
            f"|log2FC| ≥ {log2fc_cutoff})"
        ),
        xaxis_title="log2 Fold Change",
        yaxis_title=(
            f"-log10({significance_column})"
        ),
        legend_title=""
    )

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
    """
    Create gene expression boxplot(s).
    """

    gene_col = expression_df.columns[0]

    expr = expression_df.set_index(
        gene_col
    )

    expr = expr.apply(
        pd.to_numeric,
        errors="coerce"
    )

    if apply_log2:
        expr = np.log2(expr + 1)

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
        on="Sample",
        how="left"
    )
    
    groups = plot_df[group_column].dropna().unique()

    color_map = {
        group: TAB10[i % len(TAB10)]
        for i, group in enumerate(groups)
    }
    
    if plot_mode == "Single Gene":

        plot_df = plot_df[
            plot_df["Gene"] == selected_gene
        ]

        fig = px.box(
            plot_df,
            x=group_column,
            y="Expression",
            color=group_column,
            points="all",
            title=f"{selected_gene} Expression",
            color_discrete_map=color_map
        )

    else:

        fig = px.box(
            plot_df,
            x="Gene",
            y="Expression",
            color=group_column,
            points="all",
            color_discrete_map=color_map
        )



    #
    # Center points on boxplots
    #
    fig.update_traces(
        pointpos=0,
        jitter=0.08,
        marker=dict(
            size=3,
        ),
        selector=dict(
            type="box"
        )
    )

    fig.update_layout(
        template="plotly_white",
        width=width,
        height=height
    )

    if y_range is not None:

        fig.update_yaxes(
            range=y_range
        )

    return fig