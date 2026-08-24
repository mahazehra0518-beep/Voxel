"""
Unit tests for EpitopeDatasetEngine, 3-tier filtering, exclusion, rejection tracking,
and prioritization in VOXEL.
"""

import io
import pandas as pd
import pytest
from core.epitope_engine import EpitopeDatasetEngine, EpitopeFilteringConfig


@pytest.fixture
def sample_epitope_df():
    data = {
        "Peptide": [
            "YLQPRTFLL",   # Strong binder, safe, antigenic -> PASS
            "CYGVSPTKL",   # Strong binder, but TOXIC -> FAIL (Toxic)
            "GVSPTKLND",   # Good binder, but ALLERGENIC -> FAIL (Allergen)
            "WNSNNLDSK",   # Weak binder (IC50 = 1500) -> FAIL (IC50 & Rank)
            "KIADYNYKL",   # Strong binder, safe, antigenic -> PASS
            "LDSKVGGNY",   # Very weak binder (IC50 = 3500) -> FAIL (IC50)
            "NATRFASVY"    # Strong binder, safe -> PASS
        ],
        "MHC_Allele": [
            "HLA-A*02:01",
            "HLA-A*02:01",
            "HLA-A*02:01",
            "HLA-A*01:01",
            "HLA-A*02:01",
            "HLA-A*01:01",
            "HLA-A*24:02"
        ],
        "IC50_nM": [12.4, 8.5, 45.0, 1500.0, 24.6, 3500.0, 28.7],
        "Rank_%": [0.08, 0.05, 0.30, 4.50, 0.15, 8.20, 0.18],
        "Antigenicity": [0.78, 0.71, 0.68, 0.42, 0.65, 0.31, 0.69],
        "Toxicity": ["Non-Toxic", "Toxic", "Non-Toxic", "Non-Toxic", "Non-Toxic", "Non-Toxic", "Non-Toxic"],
        "Allergenicity": ["Non-Allergen", "Non-Allergen", "Allergen", "Non-Allergen", "Non-Allergen", "Non-Allergen", "Non-Allergen"]
    }
    return pd.DataFrame(data)


def test_auto_detect_columns(sample_epitope_df):
    mapping = EpitopeDatasetEngine.auto_detect_columns(sample_epitope_df)
    assert mapping["epitope"] == "Peptide"
    assert mapping["allele"] == "MHC_Allele"
    assert mapping["ic50"] == "IC50_nM"
    assert mapping["rank"] == "Rank_%"
    assert mapping["antigenicity"] == "Antigenicity"
    assert mapping["toxicity"] == "Toxicity"
    assert mapping["allergenicity"] == "Allergenicity"


def test_filtering_and_exclusion(sample_epitope_df):
    mapping = EpitopeDatasetEngine.auto_detect_columns(sample_epitope_df)
    config = EpitopeFilteringConfig(
        enable_ic50_filter=True,
        max_ic50_nm=500.0,
        enable_rank_filter=True,
        max_binding_rank=2.0,
        exclude_toxic=True,
        exclude_allergenic=True
    )
    result = EpitopeDatasetEngine.evaluate_dataset(sample_epitope_df, mapping, config)
    
    assert result.total_candidates == 7
    # Passing candidates: YLQPRTFLL, KIADYNYKL, NATRFASVY (3 candidates)
    assert result.passed_candidates_count == 3
    assert result.rejected_candidates_count == 4
    
    assert len(result.prioritized_df) == 3
    assert "Priority_Rank" in result.prioritized_df.columns
    assert "Priority_Score" in result.prioritized_df.columns
    assert result.prioritized_df.iloc[0]["Priority_Rank"] == 1


def test_rejection_reasons_logged(sample_epitope_df):
    mapping = EpitopeDatasetEngine.auto_detect_columns(sample_epitope_df)
    config = EpitopeFilteringConfig(
        enable_ic50_filter=True,
        max_ic50_nm=500.0,
        enable_rank_filter=True,
        max_binding_rank=2.0,
        exclude_toxic=True,
        exclude_allergenic=True
    )
    result = EpitopeDatasetEngine.evaluate_dataset(sample_epitope_df, mapping, config)
    audited = result.audited_full_df

    # Check Toxic row (CYGVSPTKL)
    toxic_row = audited[audited["Peptide"] == "CYGVSPTKL"].iloc[0]
    assert toxic_row["Screening_Status"] == "REJECTED"
    assert "Toxic" in toxic_row["Rejection_Reason"]

    # Check Allergenic row (GVSPTKLND)
    al_row = audited[audited["Peptide"] == "GVSPTKLND"].iloc[0]
    assert al_row["Screening_Status"] == "REJECTED"
    assert "Allergenic" in al_row["Rejection_Reason"]

    # Check High IC50 row (WNSNNLDSK)
    high_ic50_row = audited[audited["Peptide"] == "WNSNNLDSK"].iloc[0]
    assert high_ic50_row["Screening_Status"] == "REJECTED"
    assert "IC50" in high_ic50_row["Rejection_Reason"]

    # Check Passed row (YLQPRTFLL)
    pass_row = audited[audited["Peptide"] == "YLQPRTFLL"].iloc[0]
    assert pass_row["Screening_Status"] == "PASSED"
    assert "PASSED" in pass_row["Rejection_Reason"]


def test_numeric_toxicity_and_allergenicity_thresholds():
    df = pd.DataFrame({
        "Peptide": ["YLQPRTFLL", "CYGVSPTKL"],
        "MHC": ["HLA-A*02:01", "HLA-A*02:01"],
        "IC50": [20.0, 30.0],
        "Rank": [0.1, 0.2],
        "Toxin_Score": [0.15, 0.85],
        "Allergen_Score": [0.10, 0.20]
    })
    mapping = EpitopeDatasetEngine.auto_detect_columns(df)
    config = EpitopeFilteringConfig(
        exclude_toxic=True,
        toxin_score_cutoff=0.60,
        exclude_allergenic=True,
        allergen_score_cutoff=0.50
    )
    result = EpitopeDatasetEngine.evaluate_dataset(df, mapping, config)
    assert result.passed_candidates_count == 1
    assert result.rejected_candidates_count == 1
    assert "CYGVSPTKL" not in result.prioritized_df["Peptide"].values
