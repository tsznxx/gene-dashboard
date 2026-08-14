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

    fig.update_traces(
        marker=dict(
            size=10
        )
    )

    fig.update_layout(
        height=700,
        template="plotly_white",
        xaxis_title=(
            f"PC1 ({explained_variance[0]*100:.1f}%)"
        ),
        yaxis_title=(
            f"PC2 ({explained_variance[1]*100:.1f}%)"
        )
    )

    return fig