# visualization.py

import plotly.express as px


def create_pca_plot(
    pca_df,
    color_column,
    explained_variance
):

    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color=color_column,
        hover_data=["Sample"],
        title="Principal Component Analysis"
    )

    fig.update_layout(
        xaxis_title=(
            f"PC1 "
            f"({explained_variance[0]*100:.1f}%)"
        ),
        yaxis_title=(
            f"PC2 "
            f"({explained_variance[1]*100:.1f}%)"
        ),
        height=700
    )

    return fig