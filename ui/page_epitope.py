"""
Module 2: Epitope Analysis, Screening, Exclusion, and Prioritization UI for VOXEL.
"""

import os
import streamlit as st
import pandas as pd
from core.epitope_engine import EpitopeDatasetEngine, EpitopeFilteringConfig
from core.exporters import DataExporter
from ui.theme import get_badge_html, get_disclaimer_html
from ui.components import render_header, create_rejection_donut_chart, create_ic50_distribution_chart


def load_sample_dataset(dataset_type: str) -> pd.DataFrame:
    path_map = {
        "IEDB MHC-I Epitopes (CSV)": "sample_data/iedb_mhc1_predictions.csv",
        "NetMHCpan Output (TSV)": "sample_data/netmhcpan_output.tsv",
        "Comprehensive Multi-Criteria Dataset (Excel)": "sample_data/comprehensive_epitopes.xlsx"
    }
    file_path = path_map.get(dataset_type)
    if file_path and os.path.exists(file_path):
        if file_path.endswith(".csv"):
            return pd.read_csv(file_path)
        elif file_path.endswith(".tsv"):
            return pd.read_csv(file_path, sep="\t")
        elif file_path.endswith(".xlsx"):
            return pd.read_excel(file_path)
    return pd.DataFrame()


def render_epitope_page():
    render_header(
        title="MODULE 2: Epitope Screening & Prioritization",
        subtitle="Multi-criteria filtering, exclusion, audit trail, and prioritization of predicted epitopes",
        badge_text="Post-Prediction Decision Engine"
    )

    # 1. Upload / Ingestion Section
    st.markdown("### 📥 1. Upload Epitope Prediction Dataset")
    
    col_upload1, col_upload2 = st.columns([1.5, 1])
    with col_upload2:
        sample_choice = st.selectbox(
            "Or Load Quick Sample Dataset",
            ["-- Select sample --", "IEDB MHC-I Epitopes (CSV)", "NetMHCpan Output (TSV)", "Comprehensive Multi-Criteria Dataset (Excel)"]
        )

    with col_upload1:
        uploaded_file = st.file_uploader(
            "Upload prediction output file (CSV, TSV, or Excel)",
            type=["csv", "tsv", "txt", "xlsx", "xls"]
        )

    raw_df = None
    if uploaded_file is not None:
        try:
            raw_df = EpitopeDatasetEngine.parse_dataset(uploaded_file, uploaded_file.name)
            st.success(f"Successfully loaded `{uploaded_file.name}` ({len(raw_df)} rows, {len(raw_df.columns)} columns)")
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
            return
    elif sample_choice != "-- Select sample --":
        raw_df = load_sample_dataset(sample_choice)
        st.info(f"Loaded sample dataset: `{sample_choice}` ({len(raw_df)} candidates)")

    if raw_df is None or raw_df.empty:
        st.info("👆 Please upload a prediction dataset (or select a sample) above to proceed.")
        return

    # Dataset Preview & Column Mapping
    with st.expander("🔍 Dataset Preview & Column Mapping", expanded=False):
        st.dataframe(raw_df.head(5), use_container_width=True)
        st.markdown("#### Automated Column Mapping")
        st.caption("VOXEL auto-detects column names. You can adjust the mappings if needed:")
        
        auto_map = EpitopeDatasetEngine.auto_detect_columns(raw_df)
        all_cols = ["(None)"] + list(raw_df.columns)

        def get_col_index(std_key: str):
            val = auto_map.get(std_key)
            return all_cols.index(val) if val in all_cols else 0

        map_cols1, map_cols2 = st.columns(2)
        with map_cols1:
            sel_epitope = st.selectbox("Epitope / Peptide Column*", all_cols, index=get_col_index("epitope"))
            sel_allele = st.selectbox("MHC / HLA Allele Column", all_cols, index=get_col_index("allele"))
            sel_ic50 = st.selectbox("IC50 Binding Affinity Column", all_cols, index=get_col_index("ic50"))
            sel_rank = st.selectbox("Percentile Binding Rank Column", all_cols, index=get_col_index("rank"))
        with map_cols2:
            sel_ag = st.selectbox("Antigenicity Column", all_cols, index=get_col_index("antigenicity"))
            sel_tx = st.selectbox("Toxicity Column", all_cols, index=get_col_index("toxicity"))
            sel_al = st.selectbox("Allergenicity Column", all_cols, index=get_col_index("allergenicity"))

        active_mapping = {
            "epitope": None if sel_epitope == "(None)" else sel_epitope,
            "allele": None if sel_allele == "(None)" else sel_allele,
            "ic50": None if sel_ic50 == "(None)" else sel_ic50,
            "rank": None if sel_rank == "(None)" else sel_rank,
            "antigenicity": None if sel_ag == "(None)" else sel_ag,
            "toxicity": None if sel_tx == "(None)" else sel_tx,
            "allergenicity": None if sel_al == "(None)" else sel_al
        }

    if not active_mapping.get("epitope"):
        st.error("Please specify a valid Epitope/Peptide sequence column to continue.")
        return

    # 2. Multi-Criteria Configuration
    st.markdown("### ⚙️ 2. Multi-Criteria Screening & Prioritization Configuration")
    
    tab_filter, tab_exclude, tab_prioritize = st.tabs([
        "1. 📊 Filtering Criteria (Thresholds)",
        "2. 🚫 Exclusion Criteria (Safety Filters)",
        "3. ⭐ Prioritization Criteria (MCDA Weights)"
    ])

    with tab_filter:
        st.markdown("##### Continuous Quantitative Thresholds")
        f_col1, f_col2, f_col3 = st.columns(3)
        
        with f_col1:
            enable_ic50 = st.checkbox("Filter by IC50 Binding Affinity", value=True)
            preset_ic50 = st.radio(
                "Standard IC50 Preset",
                ["Strong Binder (<50 nM)", "Intermediate Binder (<500 nM)", "Weak Binder (<5000 nM)", "Custom"],
                index=1
            )
            if preset_ic50 == "Strong Binder (<50 nM)":
                max_ic50 = 50.0
            elif preset_ic50 == "Intermediate Binder (<500 nM)":
                max_ic50 = 500.0
            elif preset_ic50 == "Weak Binder (<5000 nM)":
                max_ic50 = 5000.0
            else:
                max_ic50 = st.number_input("Custom Max IC50 (nM)", min_value=1.0, max_value=50000.0, value=500.0, step=10.0)

        with f_col2:
            enable_rank = st.checkbox("Filter by Percentile Binding Rank", value=True)
            max_rank = st.slider("Max Percentile Rank Cutoff (%)", min_value=0.1, max_value=10.0, value=2.0, step=0.1,
                                 help="Lower percentile rank indicates stronger relative binding affinity.")

        with f_col3:
            enable_ag = st.checkbox("Filter by Antigenicity Score", value=False)
            min_ag = st.slider("Min Antigenicity Score", min_value=0.10, max_value=1.0, value=0.50, step=0.05)

    with tab_exclude:
        st.markdown("##### Binary / Safety Exclusion Filters")
        e_col1, e_col2, e_col3 = st.columns(3)
        
        with e_col1:
            exclude_tx = st.checkbox("Exclude Toxic Candidates", value=True, help="Removes candidates flagged as toxic or with high toxicity score.")
            tx_cutoff = st.slider("Toxin Score Cutoff", min_value=0.1, max_value=1.0, value=0.60, step=0.05)

        with e_col2:
            exclude_al = st.checkbox("Exclude Allergenic Candidates", value=True, help="Removes candidates flagged as allergenic.")
            al_cutoff = st.slider("Allergen Score Cutoff", min_value=0.1, max_value=1.0, value=0.50, step=0.05)

        with e_col3:
            exclude_amb = st.checkbox("Exclude Non-Standard / Ambiguous AAs", value=True, help="Removes peptides containing X, B, Z, J, etc.")

    with tab_prioritize:
        st.markdown("##### Multi-Criteria Decision Analysis (MCDA) Scoring Weights")
        st.caption("Relative weight assigned to each component when computing the Priority Score (0-100) for candidates that pass all filters:")
        
        w_col1, w_col2, w_col3, w_col4 = st.columns(4)
        with w_col1:
            w_ic50 = st.slider("IC50 Affinity Weight", min_value=0.0, max_value=1.0, value=0.35, step=0.05)
        with w_col2:
            w_rank = st.slider("Binding Rank Weight", min_value=0.0, max_value=1.0, value=0.25, step=0.05)
        with w_col3:
            w_ag = st.slider("Antigenicity Weight", min_value=0.0, max_value=1.0, value=0.25, step=0.05)
        with w_col4:
            w_prom = st.slider("Allele Promiscuity Weight", min_value=0.0, max_value=1.0, value=0.15, step=0.05)

    # Build Config
    config = EpitopeFilteringConfig(
        enable_ic50_filter=enable_ic50,
        max_ic50_nm=max_ic50,
        enable_rank_filter=enable_rank,
        max_binding_rank=max_rank,
        enable_antigenicity_filter=enable_ag,
        min_antigenicity_score=min_ag,
        exclude_toxic=exclude_tx,
        toxin_score_cutoff=tx_cutoff,
        exclude_allergenic=exclude_al,
        allergen_score_cutoff=al_cutoff,
        exclude_ambiguous_aa=exclude_amb,
        weight_ic50=w_ic50,
        weight_rank=w_rank,
        weight_antigenicity=w_ag,
        weight_promiscuity=w_prom
    )

    # Execute Screening Engine
    result = EpitopeDatasetEngine.evaluate_dataset(raw_df, active_mapping, config)

    # 3. Screening Results Dashboard
    st.markdown("---")
    st.markdown("### 📊 3. Screening & Prioritization Results")

    # Metrics Row
    m_col1, m_col2, m_col3, m_col4 = st.columns(4)
    m_col1.metric("Total Candidates", f"{result.total_candidates}")
    m_col2.metric("Passed All Criteria", f"{result.passed_candidates_count}", delta=f"{(result.passed_candidates_count / max(1, result.total_candidates)) * 100:.1f}%")
    m_col3.metric("Rejected Candidates", f"{result.rejected_candidates_count}")
    m_col4.metric("Screening Efficiency", f"{100.0 - (result.passed_candidates_count / max(1, result.total_candidates)) * 100:.1f}% Filtered")

    # Visualizations Row
    v_col1, v_col2 = st.columns(2)
    with v_col1:
        donut_fig = create_rejection_donut_chart(result.rejection_breakdown)
        if donut_fig:
            st.plotly_chart(donut_fig, use_container_width=True)
        else:
            st.success("🎉 All candidates passed without any rejections!")

    with v_col2:
        if active_mapping.get("ic50") and active_mapping.get("ic50") in raw_df.columns:
            hist_fig = create_ic50_distribution_chart(raw_df, active_mapping["ic50"], max_ic50)
            st.plotly_chart(hist_fig, use_container_width=True)

    # 4. Final Prioritized List
    st.markdown("### 🏆 Final Prioritized Epitope Candidates")
    st.caption("Showing candidates that successfully satisfied all filtering and exclusion criteria, ranked by composite MCDA priority score:")

    if not result.prioritized_df.empty:
        # Display key columns first
        display_cols = ["Priority_Rank", "Priority_Score", active_mapping["epitope"]]
        if active_mapping.get("allele"):
            display_cols.append(active_mapping["allele"])
        if active_mapping.get("ic50"):
            display_cols.append(active_mapping["ic50"])
        if active_mapping.get("rank"):
            display_cols.append(active_mapping["rank"])
        display_cols.append("Allele_Promiscuity_Count")
        display_cols.append("Binding_Affinity_Category")
        
        # Add remaining columns
        other_cols = [c for c in result.prioritized_df.columns if c not in display_cols and c not in ["Screening_Status", "Rejection_Reason"]]
        final_view_cols = display_cols + other_cols

        st.dataframe(
            result.prioritized_df[final_view_cols].style.background_gradient(subset=["Priority_Score"], cmap="Blues"),
            use_container_width=True,
            hide_index=True
        )
    else:
        st.warning("⚠️ No candidates passed the current criteria combination. Try relaxing the IC50 or Rank thresholds.")

    # 5. Rejection Audit Trail
    with st.expander("📋 Detailed Rejection Audit Trail (Why each candidate was rejected)", expanded=False):
        st.caption("Complete transparency log documenting the exact biological/computational reasons for exclusion:")
        audit_display_cols = [active_mapping["epitope"]]
        if active_mapping.get("allele"):
            audit_display_cols.append(active_mapping["allele"])
        audit_display_cols.extend(["Screening_Status", "Rejection_Reason"])
        if active_mapping.get("ic50"):
            audit_display_cols.append(active_mapping["ic50"])
        
        rejected_only = result.audited_full_df[result.audited_full_df["Screening_Status"] == "REJECTED"]
        st.dataframe(rejected_only[audit_display_cols], use_container_width=True, hide_index=True)

    # 6. Export Hub
    st.markdown("### 💾 4. Download Results")
    d_col1, d_col2, d_col3 = st.columns(3)

    with d_col1:
        if not result.prioritized_df.empty:
            csv_prioritized = DataExporter.to_csv(result.prioritized_df)
            st.download_button(
                label="📥 Download Prioritized Candidates (CSV)",
                data=csv_prioritized,
                file_name="voxel_prioritized_epitopes.csv",
                mime="text/csv",
                type="primary",
                use_container_width=True
            )

    with d_col2:
        csv_audited = DataExporter.to_csv(result.audited_full_df)
        st.download_button(
            label="📥 Download Full Audited Dataset (CSV)",
            data=csv_audited,
            file_name="voxel_full_audited_epitopes.csv",
            mime="text/csv",
            use_container_width=True
        )

    with d_col3:
        excel_bytes = DataExporter.to_excel_report(result)
        st.download_button(
            label="📊 Download Complete Excel Workbook (.xlsx)",
            data=excel_bytes,
            file_name="voxel_epitope_screening_report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True
        )

    # Scientific Disclaimer
    st.markdown(get_disclaimer_html(), unsafe_allow_html=True)
