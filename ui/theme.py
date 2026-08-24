"""
VOXEL Scientific Theme and Custom Styling.
Provides clean, modern bioinformatic styles, status badges, and layout aesthetics.
"""

CUSTOM_CSS = """
<style>
    /* Main container styling */
    .main-title {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        font-weight: 800;
        font-size: 2.3rem;
        background: linear-gradient(135deg, #0ea5e9 0%, #2563eb 50%, #4f46e5 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.2rem;
    }
    
    .sub-title {
        font-size: 1.05rem;
        color: #64748b;
        margin-bottom: 1.5rem;
        font-weight: 400;
    }

    /* Metric card */
    .voxel-card {
        background-color: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 1.2rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        margin-bottom: 1rem;
    }

    /* Status badge styles */
    .badge-pass {
        background-color: #dcfce7;
        color: #15803d;
        padding: 0.25rem 0.65rem;
        border-radius: 9999px;
        font-weight: 600;
        font-size: 0.85rem;
        display: inline-block;
    }

    .badge-fail {
        background-color: #fee2e2;
        color: #b91c1c;
        padding: 0.25rem 0.65rem;
        border-radius: 9999px;
        font-weight: 600;
        font-size: 0.85rem;
        display: inline-block;
    }

    .badge-warning {
        background-color: #fef3c7;
        color: #b45309;
        padding: 0.25rem 0.65rem;
        border-radius: 9999px;
        font-weight: 600;
        font-size: 0.85rem;
        display: inline-block;
    }

    .badge-info {
        background-color: #e0f2fe;
        color: #0369a1;
        padding: 0.25rem 0.65rem;
        border-radius: 9999px;
        font-weight: 600;
        font-size: 0.85rem;
        display: inline-block;
    }

    /* Science alert callout */
    .science-disclaimer {
        background-color: #f8fafc;
        border-left: 4px solid #3b82f6;
        padding: 0.85rem 1.1rem;
        border-radius: 0 8px 8px 0;
        margin-top: 1rem;
        margin-bottom: 1rem;
        font-size: 0.9rem;
        color: #334155;
    }
</style>
"""


def get_badge_html(text: str, badge_type: str = "pass") -> str:
    """Generates styled HTML badge."""
    css_cls = f"badge-{badge_type}"
    return f'<span class="{css_cls}">{text}</span>'


def get_disclaimer_html(custom_msg: str = "") -> str:
    """Generates standard scientific disclaimer box."""
    msg = custom_msg or (
        "<strong>Notice for Vaccine Researchers:</strong> In silico predictions are computational "
        "decision-support estimates. Biological activity, safety, and immunogenicity require confirmation "
        "via laboratory assays (e.g., HLA-binding ELISA, ELISpot, cell culture, in vivo models)."
    )
    return f'<div class="science-disclaimer">{msg}</div>'
