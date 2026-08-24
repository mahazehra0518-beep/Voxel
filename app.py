"""
VOXEL — Web-Based Protein and Epitope Analysis Platform.
Main Streamlit Application Entrypoint.
"""

import streamlit as st
from ui.theme import CUSTOM_CSS
from ui.page_protein import render_protein_page
from ui.page_epitope import render_epitope_page
from ui.page_documentation import render_documentation_page

# Page Configuration
st.set_page_config(
    page_title="VOXEL | Protein & Epitope Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Apply Scientific Theme CSS
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


def main():
    # Sidebar Navigation
    st.sidebar.markdown("""
        <div style="text-align: center; padding: 10px 0;">
            <h1 style="font-size: 2.2rem; margin: 0; color: #0284c7;">🧬 VOXEL</h1>
            <p style="font-size: 0.85rem; color: #64748b; margin-top: 4px;">Vaccine Decision-Support Platform</p>
        </div>
    """, unsafe_allow_html=True)

    st.sidebar.markdown("---")
    
    app_mode = st.sidebar.radio(
        "Navigation",
        [
            "🧬 Module 1: Protein Analysis",
            "🎯 Module 2: Epitope Analysis",
            "📖 Methodology & Documentation"
        ],
        index=0
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("""
        <div style="font-size: 0.82rem; color: #64748b;">
            <b>VOXEL v1.0.0</b><br>
            Computational decision-support system for reverse vaccinology & immunoinformatics.<br><br>
            <i>Notice: For research use only. In silico estimates require laboratory validation.</i>
        </div>
    """, unsafe_allow_html=True)

    # Routing
    if "Module 1" in app_mode:
        render_protein_page()
    elif "Module 2" in app_mode:
        render_epitope_page()
    elif "Documentation" in app_mode:
        render_documentation_page()


if __name__ == "__main__":
    main()
