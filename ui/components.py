"""
Reusable UI Components and Plotly Visualizations for VOXEL.
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from typing import Dict, Any, Optional
from core.sequence_parser import PhysicochemicalProperties


def render_header(title: str, subtitle: str, badge_text: Optional[str] = "Decision Support Platform"):
    badge_html = f'<span class="badge-info">{badge_text}</span>' if badge_text else ""
    st.markdown(f"""
        <div>
            <div style="display: flex; align-items: center; gap: 10px;">
                <span class="main-title">{title}</span>
                {badge_html}
            </div>
            <p class="sub-title">{subtitle}</p>
        </div>
    """, unsafe_allow_html=True)


def create_radar_chart(properties: PhysicochemicalProperties) -> go.Figure:
    """Creates a radar chart summarizing key sequence attributes."""
    categories = [
        "Helix Fraction (%)",
        "Turn Fraction (%)",
        "Sheet Fraction (%)",
        "Aromaticity (x100)",
        "pI (x10)"
    ]
    values = [
        properties.secondary_structure_fractions.get("Helix", 0),
        properties.secondary_structure_fractions.get("Turn", 0),
        properties.secondary_structure_fractions.get("Sheet", 0),
        properties.aromaticity * 100,
        properties.isoelectric_point * 10
    ]
    # Close polygon
    categories.append(categories[0])
    values.append(values[0])

    fig = go.Figure(data=go.Scatterpolar(
        r=values,
        theta=categories,
        fill="toself",
        fillcolor="rgba(14, 165, 233, 0.25)",
        line=dict(color="#0284c7", width=2)
    ))
    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 100], color="#64748b")
        ),
        margin=dict(l=40, r=40, t=30, b=30),
        height=320,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )
    return fig


def create_rejection_donut_chart(breakdown: Dict[str, int]) -> Optional[go.Figure]:
    """Generates an interactive donut chart showing the breakdown of rejection causes."""
    filtered_data = {k: v for k, v in breakdown.items() if v > 0}
    if not filtered_data:
        return None

    labels = list(filtered_data.keys())
    values = list(filtered_data.values())

    colors = ["#ef4444", "#f97316", "#eab308", "#8b5cf6", "#ec4899", "#64748b"]

    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.55,
        marker=dict(colors=colors[:len(labels)]),
        textinfo="label+value+percent",
        insidetextorientation="radial"
    )])
    fig.update_layout(
        title="<b>Rejection Bottleneck Causes</b>",
        margin=dict(l=20, r=20, t=40, b=20),
        height=340,
        showlegend=False,
        paper_bgcolor="rgba(0,0,0,0)"
    )
    return fig


def create_ic50_distribution_chart(df: pd.DataFrame, ic50_col: str, max_cutoff: float) -> go.Figure:
    """Creates a log-scale histogram of candidate IC50 affinities with threshold indicator."""
    valid_ic50 = pd.to_numeric(df[ic50_col], errors="coerce").dropna()
    valid_ic50 = valid_ic50[valid_ic50 > 0]

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=valid_ic50,
        nbinsx=30,
        marker_color="#0284c7",
        opacity=0.75,
        name="Candidates"
    ))
    # Threshold vertical line
    fig.add_vline(
        x=max_cutoff,
        line_width=2.5,
        line_dash="dash",
        line_color="#dc2626",
        annotation_text=f"Cutoff: {max_cutoff} nM",
        annotation_position="top right"
    )
    fig.update_layout(
        title="<b>IC50 Affinity Distribution (nM)</b>",
        xaxis_title="IC50 (nM)",
        yaxis_title="Count",
        xaxis_type="log",
        margin=dict(l=40, r=40, t=40, b=40),
        height=320,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )
    return fig
