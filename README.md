# gene-dashboard

A Streamlit app to upload gene expression tables and run basic analyses and visualizations (volcano plot, heatmap, survival plot, PCA, boxplot).

This repo initializes a simple, extendable Streamlit dashboard.

Usage

1. Install dependencies (create a virtualenv first):

```bash
pip install -r requirements.txt
```

2. Run the app:

```bash
streamlit run app.py
```

3. Upload a gene expression table (CSV or TSV). Expected formats (examples):

- Genes as rows: first column contains gene names and the remaining columns are samples with expression values.
- Genes as columns: check the "Transpose data" option in the app.

Optional: upload a sample metadata file (CSV/TSV) with sample IDs, group labels, or survival columns (time, event) to enable group comparisons and Kaplan–Meier plots.

What's included

- app.py: Streamlit application with upload and basic analysis tools
- requirements.txt: necessary Python packages
- .gitignore: standard Python ignores

Extend

This is a starting point—add more visualizations, advanced stats, and improved input validation as needed.
