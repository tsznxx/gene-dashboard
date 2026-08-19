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


def apply_publication_style(
    fig,
    width
):

    legend_gap_px = 40

    legend_x = (
        1
        + legend_gap_px / width
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
    fig = apply_publication_style(fig,width)
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
    fig = apply_publication_style(fig,width)
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
    fig = apply_publication_style(fig,width)
    return fig