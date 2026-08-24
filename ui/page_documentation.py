"""
Documentation, Methodology, and Threshold Reference Guide for VOXEL.
"""

import streamlit as st
from ui.components import render_header
from ui.theme import get_disclaimer_html


def render_documentation_page():
    render_header(
        title="VOXEL Scientific Documentation & Methodology Guide",
        subtitle="Computational decision-support framework for vaccine antigen & epitope design",
        badge_text="Scientific Reference"
    )

    st.markdown("""
    ### 🔬 1. Platform Overview & Workflow

    **VOXEL** bridges two critical phases in reverse vaccinology and immunoinformatics pipelines:

    1. **Pre-Epitope Protein Screening (Module 1)**:  
       Before running computationally intensive epitope binding algorithms, whole candidate protein antigens must be screened to verify:
       - **Protective Antigenicity**: Likelihood of inducing humoral/cellular protective immunity.
       - **Safety & Non-Allergenicity**: Freedom from cross-reactive IgE hypersensitivity epitopes.
       - **Non-Toxicity**: Absence of cytotoxic, hemolytic, or neurotoxic motifs.

    2. **Post-Prediction Epitope Decision Engine (Module 2)**:  
       MHC prediction tools (e.g. IEDB, NetMHCpan) generate hundreds of candidate peptides. VOXEL applies a **3-Tier Decision Support System**:
       - **Tier 1 — Continuous Filtering Criteria**: IC50 binding affinity ranges, percentile binding rank, minimum antigenicity.
       - **Tier 2 — Safety Exclusion Criteria**: Strict elimination of toxic, allergenic, or ambiguous sequences.
       - **Tier 3 — Multi-Criteria Prioritization (MCDA)**: Normalized ranking accounting for binding affinity, rank, antigenicity, and multi-allele promiscuity.
       - **Rejection Audit Trail**: Every rejected candidate is explicitly documented with its exact failure causes.

    ---

    ### 📐 2. Biological Threshold Reference Standards

    | Metric | Standard Default Threshold | Biological Interpretation | Reference / Standard |
    | :--- | :--- | :--- | :--- |
    | **$IC_{50}$ (Strong Binder)** | $< 50\text{ nM}$ | High-affinity binding to MHC groove | Sette et al. / IEDB Standard |
    | **$IC_{50}$ (Intermediate)** | $50 - 500\text{ nM}$ | Threshold generally associated with immunogenicity | Sette et al. (1994) |
    | **$IC_{50}$ (Weak Binder)** | $500 - 5000\text{ nM}$ | Borderline binding affinity | IEDB Analysis Resource |
    | **Binding Rank (%Rank)** | $\le 0.5\%$ (Strong) / $\le 2.0\%$ (Weak) | Percentile rank against random natural peptides | NetMHCpan 4.1 / Nielsen et al. |
    | **VaxiJen (Bacteria / Virus)** | Score $\ge 0.40$ | Protective antigen classification | Doytchinova & Flower (2007) |
    | **VaxiJen (Tumour / Fungi)** | Score $\ge 0.50$ | Protective / tumor antigen threshold | Doytchinova & Flower (2007) |
    | **ToxinPred2 Cutoff** | Score $\ge 0.60$ | High probability of peptide toxicity | Sharma et al. (2022) |
    | **FAO/WHO Allergen Rule** | 6-mer exact match / >35% identity | Potential for IgE cross-reactivity | FAO/WHO Codex Alimentarius (2001) |

    ---

    ### 📚 3. Scientific Citations & Model Attribution

    1. **VaxiJen v2.0**:
       > Doytchinova, I. A., & Flower, D. R. (2007). *VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines.* **BMC Bioinformatics**, 8(1), 4.
    
    2. **ToxinPred2**:
       > Sharma, N., Naorem, L. D., Jain, S., & Raghava, G. P. (2022). *ToxinPred2: an improved method for predicting toxicity of proteins.* **Briefings in Bioinformatics**, 23(5), bbac174.

    3. **AllerTOP v2.0 / v2.1**:
       > Dimitrov, I., Bangov, I., Flower, D. R., & Doytchinova, I. (2014). *AllerTOP v.2—a server for in silico prediction of allergens.* **Journal of Molecular Modeling**, 20(6), 2278.

    4. **Kolaskar & Tongaonkar Antigenic Index**:
       > Kolaskar, A. S., & Tongaonkar, P. C. (1990). *A semi-empirical method for prediction of antigenic determinants on protein antigens.* **FEBS Letters**, 276(1-2), 172-174.

    5. **FAO/WHO Allergenicity Bioinformatic Guideline**:
       > Food and Agriculture Organization / World Health Organization. (2001). *Evaluation of Allergenicity of Genetically Modified Foods.* Report of a Joint FAO/WHO Expert Consultation.
    """)

    # Scientific Disclaimer Box
    st.markdown(get_disclaimer_html(), unsafe_allow_html=True)
