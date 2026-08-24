"""
Unit tests for ProteinAnalyzer and biological model adapters in VOXEL.
"""

import pytest
from core.sequence_parser import SequenceParser
from core.protein_analyzer import ProteinAnalyzer
from models.base import PredictionCategory, PredictorStatus
from models.antigenicity.vaxijen_adapter import VaxiJenAdapter
from models.toxicity.toxinpred_adapter import ToxinPredAdapter
from models.allergenicity.allertop_adapter import AllerTOPAdapter


def test_vaxijen_adapter_user_score():
    adapter = VaxiJenAdapter()
    res = adapter.predict("ACDEFGHIKLMNPQRSTVWY", user_score=0.65, organism="Virus")
    assert res.category == PredictionCategory.ANTIGENICITY
    assert res.prediction_score == 0.65
    assert res.classification == "Probable Antigen"
    assert res.is_favorable_for_vaccine() is True
    assert "0.65" in res.interpretation
    assert "experimental validation" in res.disclaimer.lower()


def test_toxinpred_adapter_user_score():
    adapter = ToxinPredAdapter()
    res_toxic = adapter.predict("ACDEFGHIKLMNPQRSTVWY", user_score=0.85, custom_threshold=0.60)
    assert res_toxic.category == PredictionCategory.TOXICITY
    assert res_toxic.prediction_score == 0.85
    assert res_toxic.classification == "Toxic"
    assert res_toxic.is_favorable_for_vaccine() is False

    res_safe = adapter.predict("ACDEFGHIKLMNPQRSTVWY", user_score=0.20, custom_threshold=0.60)
    assert res_safe.classification == "Non-Toxic"
    assert res_safe.is_favorable_for_vaccine() is True


def test_allertop_adapter_user_classification():
    adapter = AllerTOPAdapter()
    res_al = adapter.predict("ACDEFGHIKLMNPQRSTVWY", user_classification="Probable Allergen")
    assert res_al.category == PredictionCategory.ALLERGENICITY
    assert res_al.classification == "Probable Allergen"
    assert res_al.is_favorable_for_vaccine() is False

    res_non = adapter.predict("ACDEFGHIKLMNPQRSTVWY", user_classification="Probable Non-Allergen")
    assert res_non.classification == "Probable Non-Allergen"
    assert res_non.is_favorable_for_vaccine() is True


def test_protein_analyzer_aggregation():
    analyzer = ProteinAnalyzer()
    records = SequenceParser.parse_fasta(">test\nRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNS\n")
    report = analyzer.analyze(
        records[0],
        organism="Virus",
        user_antigenicity_score=0.75,
        user_allergenicity_cls="Non-Allergen",
        user_toxicity_score=0.10
    )
    assert report.passes_all_criteria is True
    assert "HIGH POTENTIAL" in report.overall_suitability_verdict
    assert report.overall_suitability_score > 70.0
    assert len(report.recommendations) > 0
