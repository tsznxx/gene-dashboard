import streamlit as st

from data_loader import (
    load_expression_file,
    load_metadata_file,
    validate_expression_matrix,
    validate_metadata,
    validate_sample_matching,
    summarize_expression,
    summarize_metadata
)

st.set_page_config(
    page_title="Gene Expression Dashboard",
    layout="wide"
)

st.title("Gene Expression Dashboard")

st.header("Module 2: Data Upload & Validation")

left, right = st.columns(2)

with left:

    uploaded_expression = st.file_uploader(
        "Expression Matrix",
        type=["csv"]
    )

with right:

    uploaded_metadata = st.file_uploader(
        "Metadata Table",
        type=["csv"]
    )

if uploaded_expression and uploaded_metadata:

    try:

        expr_df = load_expression_file(
            uploaded_expression
        )

        meta_df = load_metadata_file(
            uploaded_metadata
        )

        expr_errors = validate_expression_matrix(
            expr_df
        )

        meta_errors = validate_metadata(
            meta_df
        )

        all_errors = expr_errors + meta_errors

        if all_errors:

            st.error(
                "Validation failed."
            )

            for err in all_errors:
                st.write("•", err)

            st.stop()

        match_result = validate_sample_matching(
            expr_df,
            meta_df
        )

        if not match_result["matching"]:

            st.error(
                "Sample IDs do not match."
            )

            if match_result["missing_in_metadata"]:
                st.write(
                    "Missing in metadata:"
                )
                st.write(
                    match_result["missing_in_metadata"]
                )

            if match_result["missing_in_expression"]:
                st.write(
                    "Missing in expression matrix:"
                )
                st.write(
                    match_result["missing_in_expression"]
                )

            st.stop()

        st.success(
            "Expression matrix and metadata validated."
        )

        st.session_state["expression_df"] = expr_df
        st.session_state["metadata_df"] = meta_df

        st.subheader("Dataset Summary")

        e_summary = summarize_expression(
            expr_df
        )

        m_summary = summarize_metadata(
            meta_df
        )

        c1, c2, c3 = st.columns(3)

        with c1:
            st.metric(
                "Genes",
                e_summary["Genes"]
            )

        with c2:
            st.metric(
                "Samples",
                e_summary["Samples"]
            )

        with c3:
            st.metric(
                "Groups",
                m_summary.get("Groups", "NA")
            )

        tab1, tab2 = st.tabs(
            [
                "Expression Matrix",
                "Metadata"
            ]
        )

        with tab1:
            st.dataframe(
                expr_df.head(20),
                use_container_width=True
            )

        with tab2:
            st.dataframe(
                meta_df.head(20),
                use_container_width=True
            )

    except Exception as e:

        st.exception(e)

else:

    st.info(
        "Please upload both files."
    )