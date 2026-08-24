"""
Module 1: Protein Sequence Analysis UI for VOXEL.
"""

import os
import streamlit as st
import pandas as pd
from core.sequence_parser import SequenceParser, ProteinSequence
from core.protein_analyzer import ProteinAnalyzer
from ui.theme import get_badge_html, get_disclaimer_html
from ui.components import render_header, create_radar_chart


def load_sample_sequence(sample_type: str) -> str:
    path_map = {
        "SARS-CoV-2 Spike RBD": "sample_data/spike_glycoprotein.fasta",
        "M. tuberculosis Ag85B": "sample_data/ag85b_mycobacterium.fasta"
    }
    file_path = path_map.get(sample_type)
    if file_path and os.path.exists(file_path):
        with open(file_path, "r") as f:
            return f.read()
    return ""


def render_protein_page():
    render_header(
        title="MODULE 1: Protein Sequence Analysis",
        subtitle="Pre-epitope computational screening for Antigenicity, Allergenicity, and Toxicity",
        badge_text="Pre-Screening Pipeline"
    )

    col_main, col_sidebar = st.columns([2, 1])

    with col_sidebar:
        st.markdown("### ⚙️ Analysis Parameters")
        organism = st.selectbox(
            "Target Organism Class",
            ["Virus", "Bacteria", "Tumour", "Parasite", "Fungi"],
            index=0,
            help="Selects the organism-specific model and threshold for VaxiJen."
        )

        with st.expander("🛠️ Prediction Cutoff Thresholds", expanded=False):
            ag_thresh = st.slider(
                "VaxiJen Antigenicity Cutoff",
                min_value=0.20, max_value=0.80,
                value=0.40 if organism in ["Virus", "Bacteria"] else 0.50,
                step=0.05,
                help="Scores at or above this threshold classify the sequence as a protective antigen."
            )
            tx_thresh = st.slider(
                "ToxinPred2 Toxicity Cutoff",
                min_value=0.20, max_value=0.90,
                value=0.60,
                step=0.05,
                help="Scores at or above this threshold classify the sequence as toxic."
            )
            al_thresh = st.slider(
                "Allergenicity Cutoff",
                min_value=0.20, max_value=0.90,
                value=0.50,
                step=0.05,
                help="Scores at or above this threshold classify the sequence as an allergen."
            )

        with st.expander("📝 External Server Verified Scores (Optional)", expanded=False):
            st.caption("If you ran the sequence on external official web servers, enter the verified outputs below:")
            user_ag = st.number_input("Verified VaxiJen Score", min_value=0.0, max_value=2.0, value=None, step=0.01)
            user_tx = st.number_input("Verified ToxinPred2 Score", min_value=0.0, max_value=1.0, value=None, step=0.01)
            user_al = st.selectbox("Verified AllerTOP Verdict", ["(None / Use Engine)", "Probable Non-Allergen", "Probable Allergen"], index=0)
            user_al_val = None if user_al == "(None / Use Engine)" else user_al

    with col_main:
        st.markdown("### 📥 Protein Input")
        sample_choice = st.selectbox(
            "Load Quick Sample Sequence (Optional)",
            ["-- Select a sample --", "SARS-CoV-2 Spike RBD", "M. tuberculosis Ag85B"]
        )
        
        default_text = ""
        if sample_choice != "-- Select a sample --":
            default_text = load_sample_sequence(sample_choice)

        uploaded_fasta = st.file_uploader("Upload FASTA file (.fasta, .fa, .txt)", type=["fasta", "fa", "txt"])
        if uploaded_fasta is not None:
            default_text = uploaded_fasta.read().decode("utf-8")

        fasta_input = st.text_area(
            "Or Paste FASTA / Plain Sequence Text",
            value=default_text,
            height=160,
            placeholder=">Protein_ID [Organism]\nACDEFGHIKLMNPQRSTVWY..."
        )

        analyze_btn = st.button("🚀 Run Comprehensive Protein Analysis", type="primary", use_container_width=True)

    if analyze_btn or (fasta_input and len(fasta_input.strip()) > 0):
        if not fasta_input.strip():
            st.warning("Please provide a protein sequence to begin analysis.")
            return

        # Parse FASTA
        records = SequenceParser.parse_fasta(fasta_input, default_organism=organism)
        if not records:
            st.error("Could not parse any valid sequence from input.")
            return

        protein_seq = records[0]

        if protein_seq.validation_warnings:
            for w in protein_seq.validation_warnings:
                st.warning(f"⚠️ **Sequence Notice**: {w}")

        if not protein_seq.is_valid:
            st.error("Sequence validation failed. Please review sequence content.")
            return

        # Run Analysis
        analyzer = ProteinAnalyzer()
        report = analyzer.analyze(
            protein_seq=protein_seq,
            organism=organism,
            user_antigenicity_score=user_ag,
            user_allergenicity_cls=user_al_val,
            user_toxicity_score=user_tx,
            antigenicity_threshold=ag_thresh,
            allergenicity_threshold=al_thresh,
            toxicity_threshold=tx_thresh
        )

        st.markdown("---")
        st.markdown(f"## 📊 Analysis Report: `{protein_seq.identifier}`")
        if protein_seq.organism:
            st.caption(f"**Organism Source**: {protein_seq.organism}")

        # Overall Verdict Callout
        verdict_badge = "pass" if report.passes_all_criteria else ("fail" if "TOXIC" in report.overall_suitability_verdict or "RISK" in report.overall_suitability_verdict else "warning")
        
        st.markdown(f"""
            <div class="voxel-card" style="border-left: 6px solid {'#22c55e' if report.passes_all_criteria else '#ef4444'};">
                <div style="display:flex; justify-content:space-between; align-items:center;">
                    <div>
                        <span style="font-size: 1.2rem; font-weight:700;">{report.overall_suitability_verdict}</span>
                        <div style="margin-top: 5px;">{get_badge_html(f"Suitability Score: {report.overall_suitability_score}/100", verdict_badge)}</div>
                    </div>
                </div>
            </div>
        """, unsafe_allow_html=True)

        for rec in report.recommendations:
            st.info(f"💡 **Recommendation**: {rec}")

        # Tri-Factor Assessment Cards
        st.markdown("### 🧪 Tri-Factor Biological Evaluation")
        card_cols = st.columns(3)

        # 1. Antigenicity Card
        with card_cols[0]:
            ag_res = report.antigenicity_result
            ag_badge = "pass" if ag_res.is_favorable_for_vaccine() else "fail"
            st.markdown(f"""
                <div class="voxel-card">
                    <h4>🛡️ Antigenicity</h4>
                    <div style="margin-bottom:8px;">{get_badge_html(ag_res.classification, ag_badge)}</div>
                    <p><b>Prediction Score:</b> <code>{ag_res.prediction_score if ag_res.prediction_score is not None else 'N/A'}</code></p>
                    <p><b>Threshold:</b> ≥ {ag_res.threshold_used:.2f}</p>
                    <p style="font-size:0.85rem; color:#64748b;"><b>Tool:</b> {ag_res.tool_name}</p>
                    <p style="font-size:0.85rem;">{ag_res.interpretation}</p>
                </div>
            """, unsafe_allow_html=True)

        # 2. Allergenicity Card
        with card_cols[1]:
            al_res = report.allergenicity_result
            al_badge = "pass" if al_res.is_favorable_for_vaccine() else "fail"
            st.markdown(f"""
                <div class="voxel-card">
                    <h4>🌿 Allergenicity</h4>
                    <div style="margin-bottom:8px;">{get_badge_html(al_res.classification, al_badge)}</div>
                    <p><b>Prediction Score:</b> <code>{al_res.prediction_score if al_res.prediction_score is not None else 'N/A'}</code></p>
                    <p><b>Threshold:</b> &lt; {al_res.threshold_used:.2f}</p>
                    <p style="font-size:0.85rem; color:#64748b;"><b>Tool:</b> {al_res.tool_name}</p>
                    <p style="font-size:0.85rem;">{al_res.interpretation}</p>
                </div>
            """, unsafe_allow_html=True)

        # 3. Toxicity Card
        with card_cols[2]:
            tx_res = report.toxicity_result
            tx_badge = "pass" if tx_res.is_favorable_for_vaccine() else "fail"
            st.markdown(f"""
                <div class="voxel-card">
                    <h4>☠️ Toxicity</h4>
                    <div style="margin-bottom:8px;">{get_badge_html(tx_res.classification, tx_badge)}</div>
                    <p><b>Prediction Score:</b> <code>{tx_res.prediction_score if tx_res.prediction_score is not None else 'N/A'}</code></p>
                    <p><b>Threshold:</b> &lt; {tx_res.threshold_used:.2f}</p>
                    <p style="font-size:0.85rem; color:#64748b;"><b>Tool:</b> {tx_res.tool_name}</p>
                    <p style="font-size:0.85rem;">{tx_res.interpretation}</p>
                </div>
            """, unsafe_allow_html=True)

        # Physicochemical Properties Section
        if protein_seq.properties:
            props = protein_seq.properties
            st.markdown("### 🔬 Physicochemical Characterization")
            
            p_col1, p_col2 = st.columns([1.8, 1.2])

            with p_col1:
                metric_row1 = st.columns(3)
                metric_row1[0].metric("Length", f"{props.length} AA")
                metric_row1[1].metric("Mol Weight", f"{props.molecular_weight_kda:.2f} kDa")
                metric_row1[2].metric("Isoelectric Pt (pI)", f"{props.isoelectric_point:.2f}")

                metric_row2 = st.columns(3)
                metric_row2[0].metric("GRAVY Hydropathy", f"{props.gravy_hydropathy:.3f}")
                metric_row2[1].metric("Instability Index", f"{props.instability_index:.2f}", delta="Stable in vitro" if props.is_stable_in_vitro else "Unstable")
                metric_row2[2].metric("Net Charge (pH 7)", f"{props.charge_at_ph7:.2f}")

                with st.expander("🔍 Detailed Amino Acid Frequency Table"):
                    aa_df = pd.DataFrame([
                        {"Amino Acid": aa, "Count": props.amino_acid_counts.get(aa, 0), "Percentage (%)": props.amino_acid_percentages.get(aa, 0.0)}
                        for aa in sorted("ACDEFGHIKLMNPQRSTVWY")
                    ])
                    st.dataframe(aa_df, use_container_width=True, hide_index=True)

            with p_col2:
                st.markdown("##### Secondary Structure & Properties Radar")
                fig = create_radar_chart(props)
                st.plotly_chart(fig, use_container_width=True)

        # Explicit Scientific Disclaimer Box
        st.markdown(get_disclaimer_html(), unsafe_allow_html=True)
