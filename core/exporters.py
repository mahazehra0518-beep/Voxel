"""
Data Exporters for VOXEL Platform.
Supports CSV, Excel (multi-tab with styling), and structured text/markdown reports.
"""

import io
import pandas as pd
from typing import Dict, Any
from core.epitope_engine import EpitopeAnalysisResult


class DataExporter:
    @staticmethod
    def to_csv(df: pd.DataFrame) -> bytes:
        """Converts DataFrame to UTF-8 CSV bytes."""
        return df.to_csv(index=False).encode("utf-8")

    @staticmethod
    def to_excel_report(result: EpitopeAnalysisResult, project_name: str = "VOXEL Analysis") -> bytes:
        """
        Generates a comprehensive multi-tab Excel workbook:
        - Tab 1: Prioritized Epitope Candidates (Passing all criteria)
        - Tab 2: Full Audited Dataset (With explicit rejection reasons)
        - Tab 3: Screening Summary & Rejection Breakdown
        """
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            # Tab 1: Prioritized candidates
            result.prioritized_df.to_excel(writer, sheet_name="Prioritized Candidates", index=False)
            
            # Tab 2: Full Audited Dataset
            result.audited_full_df.to_excel(writer, sheet_name="Full Audited Dataset", index=False)
            
            # Tab 3: Summary Metadata
            summary_data = [
                {"Metric": "Project Name", "Value": project_name},
                {"Metric": "Total Epitopes Evaluated", "Value": result.total_candidates},
                {"Metric": "Candidates Passed All Criteria", "Value": result.passed_candidates_count},
                {"Metric": "Candidates Rejected", "Value": result.rejected_candidates_count},
                {"Metric": "Pass Rate (%)", "Value": f"{(result.passed_candidates_count / max(1, result.total_candidates)) * 100:.1f}%"},
                {"Metric": "--- REJECTION BREAKDOWN ---", "Value": "---"},
            ]
            for reason, count in result.rejection_breakdown.items():
                summary_data.append({"Metric": f"Rejected by: {reason}", "Value": count})
            
            # Add configuration parameters
            summary_data.extend([
                {"Metric": "--- SCREENING CONFIGURATION ---", "Value": "---"},
                {"Metric": "IC50 Filter Enabled", "Value": str(result.config.enable_ic50_filter)},
                {"Metric": "Max IC50 (nM)", "Value": str(result.config.max_ic50_nm)},
                {"Metric": "Rank Filter Enabled", "Value": str(result.config.enable_rank_filter)},
                {"Metric": "Max %Rank Cutoff", "Value": str(result.config.max_binding_rank)},
                {"Metric": "Exclude Toxic Candidates", "Value": str(result.config.exclude_toxic)},
                {"Metric": "Exclude Allergenic Candidates", "Value": str(result.config.exclude_allergenic)},
                {"Metric": "Notice", "Value": "All predictions are computational in silico estimates requiring experimental validation."}
            ])
            
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name="Screening Summary", index=False)

        return output.getvalue()

    @staticmethod
    def generate_markdown_summary(result: EpitopeAnalysisResult) -> str:
        """Generates a structured markdown audit report."""
        lines = [
            "# VOXEL Epitope Screening Summary Report",
            "",
            "## Executive Summary",
            f"- **Total Candidates Analyzed**: {result.total_candidates}",
            f"- **Passed All Criteria**: {result.passed_candidates_count} ({result.passed_candidates_count / max(1, result.total_candidates) * 100:.1f}%)",
            f"- **Rejected Candidates**: {result.rejected_candidates_count}",
            "",
            "## Rejection Bottleneck Breakdown",
        ]
        for reason, count in result.rejection_breakdown.items():
            lines.append(f"- **{reason}**: {count} candidates failed")
            
        lines.extend([
            "",
            "## Applied Threshold Configuration",
            f"- Max $IC_{{50}}$: {result.config.max_ic50_nm} nM (Filter Enabled: {result.config.enable_ic50_filter})",
            f"- Max Percentile Binding Rank: {result.config.max_binding_rank}% (Filter Enabled: {result.config.enable_rank_filter})",
            f"- Toxicity Exclusion: {result.config.exclude_toxic}",
            f"- Allergenicity Exclusion: {result.config.exclude_allergenic}",
            "",
            "> **Scientific Disclaimer**: Predictions presented here are computed decision-support estimates. Laboratory validation (ELISpot, binding assays, animal models) is strictly required prior to clinical application."
        ])
        return "\n".join(lines)
