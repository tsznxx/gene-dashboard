import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from lifelines import KaplanMeierFitter
from scipy import stats
import plotly.express as px

st.set_page_config(page_title="Gene Dashboard", layout="wide")

st.title("Gene expression dashboard — basic analyses")

# Sidebar: file uploads and options
st.sidebar.header("Data upload & options")
expr_file = st.sidebar.file_uploader("Upload gene expression file (CSV/TSV)", type=["csv", "tsv", "txt"] )
meta_file = st.sidebar.file_uploader("(Optional) Sample metadata file (CSV/TSV)", type=["csv", "tsv", "txt"] )
transpose = st.sidebar.checkbox("Transpose data (genes in columns)", value=False)
index_col = st.sidebar.text_input("Gene ID column name (if present)", value="")

# Helper: read uploaded file
def read_table(uploaded):
    if uploaded is None:
        return None
    content = uploaded.read()
    try:
        # try CSV
        df = pd.read_csv(io.BytesIO(content))
    except Exception:
        try:
            df = pd.read_csv(io.BytesIO(content), sep="\t")
        except Exception:
            st.error("Could not read the uploaded file. Make sure it's CSV or TSV.")
            return None
    return df

expr_df = read_table(expr_file)
meta_df = read_table(meta_file)

if expr_df is None:
    st.info("Upload a gene expression table to get started.")
    st.stop()

# optionally set gene index
if index_col and index_col in expr_df.columns:
    expr_df = expr_df.set_index(index_col)

if transpose:
    expr_df = expr_df.T

st.subheader("Expression data preview")
st.dataframe(expr_df.head())

# Basic checks
st.sidebar.header("Analysis controls")
show_na = st.sidebar.checkbox("Show missing value summary", value=False)
if show_na:
    na_summary = expr_df.isna().sum()
    st.write("Missing values per column:")
    st.write(na_summary[na_summary>0])

# If metadata provided, preview
if meta_df is not None:
    st.subheader("Sample metadata preview")
    st.dataframe(meta_df.head())

# Analysis: PCA
st.header("PCA")
with st.expander("Run PCA on samples"):
    n_components = st.slider("Number of PCA components", min_value=2, max_value=5, value=2)
    # Prepare matrix: samples x genes
    if expr_df.shape[0] < expr_df.shape[1]:
        data_for_pca = expr_df.fillna(0).values
    else:
        # assume genes are rows, samples are columns -> transpose
        data_for_pca = expr_df.fillna(0).T.values
    pca = PCA(n_components=n_components)
    try:
        pcs = pca.fit_transform(data_for_pca)
        pc_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
        # attach sample names
        sample_names = (expr_df.columns if expr_df.shape[0] < expr_df.shape[1] else expr_df.index)
        pc_df["sample"] = sample_names
        if meta_df is not None:
            # try to merge on sample name column if present
            # assume first column of metadata is sample ID unless there's an explicit 'sample' column
            meta_index = None
            if "sample" in meta_df.columns:
                meta_index = "sample"
            else:
                meta_index = meta_df.columns[0]
            merged = pc_df.merge(meta_df, left_on="sample", right_on=meta_index, how="left")
        else:
            merged = pc_df

        if n_components >= 2:
            fig = px.scatter(merged, x="PC1", y="PC2", color=merged.columns[2] if merged.shape[1]>2 else None,
                             hover_name="sample")
            st.plotly_chart(fig, use_container_width=True)
            st.write("Explained variance ratio:", pca.explained_variance_ratio_)
    except Exception as e:
        st.error(f"PCA failed: {e}")

# Heatmap of top variable genes
st.header("Heatmap")
with st.expander("Top variable genes heatmap"):
    top_n = st.number_input("Top N variable genes", min_value=5, max_value=200, value=50)
    # compute variance per gene (assume genes are rows)
    if expr_df.shape[0] >= expr_df.shape[1]:
        gene_vars = expr_df.var(axis=1).sort_values(ascending=False)
        top_genes = gene_vars.head(top_n).index
        heat_df = expr_df.loc[top_genes].fillna(0)
    else:
        gene_vars = expr_df.var(axis=0).sort_values(ascending=False)
        top_genes = gene_vars.head(top_n).index
        heat_df = expr_df[top_genes].T.fillna(0)
    fig, ax = plt.subplots(figsize=(10, min(10, 0.2*len(top_genes)+3)))
    sns.heatmap(heat_df, cmap="vlag", ax=ax)
    st.pyplot(fig)

# Boxplot for selected gene(s)
st.header("Boxplot for gene(s)")
with st.expander("Boxplot / violin by group"):
    gene_choice = st.text_input("Enter gene name (or comma-separated list)")
    group_col = st.selectbox("Group column (from metadata)", options=[None] + (list(meta_df.columns) if meta_df is not None else []))
    if gene_choice:
        genes = [g.strip() for g in gene_choice.split(",")]
        # get expression series
        data = {}
        for g in genes:
            if g in expr_df.index:
                data[g] = expr_df.loc[g]
            elif g in expr_df.columns:
                data[g] = expr_df[g]
            else:
                st.warning(f"Gene {g} not found in expression table")
        if data:
            plot_df = pd.DataFrame(data)
            # align metadata
            if meta_df is not None and group_col in meta_df.columns:
                meta_index = meta_df.columns[0]
                plot_df = plot_df.T if plot_df.shape[0] != expr_df.shape[0] else plot_df
                # ensure rows are samples
                if plot_df.shape[0] == expr_df.shape[0]:
                    joined = plot_df.reset_index().rename(columns={"index":"sample"}).merge(meta_df, left_on="sample", right_on=meta_index, how="left")
                else:
                    # transpose and try merge on column names
                    joined = plot_df.T.reset_index().rename(columns={"index":"sample"}).merge(meta_df, left_on="sample", right_on=meta_index, how="left")
                plt.figure(figsize=(8,4))
                sns.boxplot(x=group_col, y=genes[0], data=joined)
                st.pyplot(plt)
            else:
                st.write(plot_df.describe())

# Volcano plot (require two groups)
st.header("Volcano plot (two-group comparison)")
with st.expander("Differential expression between two groups"):
    if meta_df is None:
        st.info("Upload sample metadata with a group column to enable volcano plot")
    else:
        group_column = st.selectbox("Select group column", options=[c for c in meta_df.columns])
        unique_groups = meta_df[group_column].dropna().unique()
        if len(unique_groups) >= 2:
            g1 = st.selectbox("Group 1", options=unique_groups, index=0)
            g2 = st.selectbox("Group 2", options=unique_groups, index=1)
            run_volcano = st.button("Run volcano")
            if run_volcano:
                # map samples to groups
                meta_index = meta_df.columns[0]
                sample_to_group = meta_df.set_index(meta_index)[group_column].to_dict()
                # get sample lists
                samples = list(expr_df.columns) if expr_df.shape[0] >= expr_df.shape[1] else list(expr_df.index)
                g1_samples = [s for s in samples if sample_to_group.get(s) == g1]
                g2_samples = [s for s in samples if sample_to_group.get(s) == g2]
                if len(g1_samples) < 2 or len(g2_samples) < 2:
                    st.error("Need at least 2 samples per group")
                else:
                    # compute t-tests per gene
                    if expr_df.shape[0] >= expr_df.shape[1]:
                        genes = expr_df.index
                        res = []
                        for gene in genes:
                            a = expr_df.loc[gene, g1_samples].dropna()
                            b = expr_df.loc[gene, g2_samples].dropna()
                            if len(a) < 2 or len(b) < 2:
                                continue
                            t, p = stats.ttest_ind(a, b, equal_var=False)
                            fc = a.mean() - b.mean()  # log fold not computed; using difference
                            res.append((gene, fc, p))
                        res_df = pd.DataFrame(res, columns=["gene","fc","pval"]) 
                        res_df["neglog10p"] = -np.log10(res_df["pval"]+1e-300)
                        fig = px.scatter(res_df, x="fc", y="neglog10p", hover_name="gene")
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.error("Expression layout not recognized for volcano calculation")
        else:
            st.info("Group column must have at least two distinct values")

# Survival plot
st.header("Survival / Kaplan–Meier")
with st.expander("Kaplan–Meier from metadata (requires 'time' and 'event' columns)"):
    if meta_df is None:
        st.info("Upload metadata with 'time' and 'event' columns to enable survival analysis")
    else:
        candidates = [c for c in meta_df.columns if c.lower() in ["time","duration","days","event","status","dead"]]
        st.write("Detected metadata columns: ", list(meta_df.columns))
        time_col = st.selectbox("Time column", options=[None]+list(meta_df.columns))
        event_col = st.selectbox("Event column", options=[None]+list(meta_df.columns))
        group_col = st.selectbox("Optional group column", options=[None]+list(meta_df.columns))
        if time_col and event_col:
            try:
                df = meta_df[[time_col, event_col, group_col]] if group_col else meta_df[[time_col, event_col]]
                df = df.dropna(subset=[time_col, event_col])
                kmf = KaplanMeierFitter()
                if group_col:
                    fig, ax = plt.subplots()
                    for name, grouped_df in df.groupby(group_col):
                        kmf.fit(grouped_df[time_col], event_observed=grouped_df[event_col], label=str(name))
                        kmf.plot_survival_function(ax=ax)
                    st.pyplot(fig)
                else:
                    kmf.fit(df[time_col], event_observed=df[event_col])
                    fig = kmf.plot_survival_function()
                    st.pyplot()
            except Exception as e:
                st.error(f"Survival plot failed: {e}")

st.sidebar.markdown("---")
st.sidebar.write("Ready. Extend the app by adding more analyses and input validation.")
