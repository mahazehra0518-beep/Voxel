"""
Module 1: Protein Sequence Analysis Pipeline for VOXEL.
Orchestrates tri-factor biological assessment:
1. Antigenicity (VaxiJen v2.0 / ACC Propensity)
2. Allergenicity (AllerTOP v2.1 / FAO-WHO Rule Engine)
3. Toxicity (ToxinPred2 / Physicochemical Risk Screen)

Provides comprehensive synthesis, suitability verdicts, and experimental validation warnings.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List
from core.sequence_parser import ProteinSequence
from models.base import PredictionResult, PredictionCategory, PredictorStatus
from models.antigenicity.vaxijen_adapter import VaxiJenAdapter
from models.toxicity.toxinpred_adapter import ToxinPredAdapter
from models.allergenicity.allertop_adapter import AllerTOPAdapter


@dataclass
class ProteinAnalysisReport:
    sequence_info: ProteinSequence
    antigenicity_result: PredictionResult
    allergenicity_result: PredictionResult
    toxicity_result: PredictionResult
    overall_suitability_verdict: str
    overall_suitability_score: float  # 0 to 100
    recommendations: List[str] = field(default_factory=list)
    
    @property
    def passes_all_criteria(self) -> bool:
        return (
            self.antigenicity_result.is_favorable_for_vaccine() and
            self.allergenicity_result.is_favorable_for_vaccine() and
            self.toxicity_result.is_favorable_for_vaccine()
        )


class ProteinAnalyzer:
    def __init__(self):
        self.vaxijen_adapter = VaxiJenAdapter()
        self.toxinpred_adapter = ToxinPredAdapter()
        self.allertop_adapter = AllerTOPAdapter()

    def analyze(
        self,
        protein_seq: ProteinSequence,
        organism: str = "Virus",
        user_antigenicity_score: Optional[float] = None,
        user_allergenicity_cls: Optional[str] = None,
        user_allergenicity_score: Optional[float] = None,
        user_toxicity_score: Optional[float] = None,
        antigenicity_threshold: Optional[float] = None,
        allergenicity_threshold: Optional[float] = None,
        toxicity_threshold: Optional[float] = None,
        attempt_remote_bridges: bool = False
    ) -> ProteinAnalysisReport:
        """
        Executes full tri-factor analysis on the protein sequence.
        """
        raw_seq = protein_seq.clean_sequence

        # 1. Antigenicity Evaluation
        ag_res = self.vaxijen_adapter.predict(
            sequence=raw_seq,
            user_score=user_antigenicity_score,
            organism=organism,
            custom_threshold=antigenicity_threshold,
            attempt_remote=attempt_remote_bridges
        )

        # 2. Allergenicity Evaluation
        al_res = self.allertop_adapter.predict(
            sequence=raw_seq,
            user_classification=user_allergenicity_cls,
            user_score=user_allergenicity_score,
            custom_threshold=allergenicity_threshold
        )

        # 3. Toxicity Evaluation
        tx_res = self.toxinpred_adapter.predict(
            sequence=raw_seq,
            user_score=user_toxicity_score,
            custom_threshold=toxicity_threshold,
            attempt_local=attempt_remote_bridges
        )

        # Calculate composite suitability score (0 - 100)
        # Higher antigenicity = higher score, lower allergenicity/toxicity = higher score
        ag_contrib = (ag_res.prediction_score if ag_res.prediction_score is not None else 0.5) * 40.0
        al_contrib = (1.0 - (al_res.prediction_score if al_res.prediction_score is not None else 0.5)) * 30.0
        tx_contrib = (1.0 - (tx_res.prediction_score if tx_res.prediction_score is not None else 0.5)) * 30.0
        
        suitability_score = round(max(0.0, min(100.0, ag_contrib + al_contrib + tx_contrib)), 1)

        # Build overall verdict and actionable recommendations
        recommendations = []
        is_ag = ag_res.is_favorable_for_vaccine()
        is_safe_al = al_res.is_favorable_for_vaccine()
        is_safe_tx = tx_res.is_favorable_for_vaccine()

        if is_ag and is_safe_al and is_safe_tx:
            verdict = "HIGH POTENTIAL — Highly Recommended for Epitope Prediction"
            recommendations.append("Protein satisfies all essential pre-screening criteria: Antigenic, Non-Allergenic, and Non-Toxic.")
            recommendations.append("Proceed to Module 2 (Epitope Screening) following MHC-I/MHC-II binding prediction.")
        elif not is_safe_tx:
            verdict = "HIGH RISK — Predicted Toxic"
            recommendations.append("Sequence flagged as potentially toxic. Identify and mutate specific toxic residues/motifs before in vivo use.")
        elif not is_safe_al:
            verdict = "MODERATE RISK — Potential Allergenicity Detected"
            recommendations.append("Sequence flagged with allergenic potential or matching cross-reactive IgE epitopes. Cross-reference with clinical allergen databases.")
        elif not is_ag:
            verdict = "LOW IMMUNOGENICITY — Low Antigenic Propensity"
            recommendations.append("Sequence has low predicted protective antigenicity. Consider exploring adjuvant conjugation or carrier protein fusion.")
        else:
            verdict = "BORDERLINE CANDIDATE — Secondary Review Advised"

        return ProteinAnalysisReport(
            sequence_info=protein_seq,
            antigenicity_result=ag_res,
            allergenicity_result=al_res,
            toxicity_result=tx_res,
            overall_suitability_verdict=verdict,
            overall_suitability_score=suitability_score,
            recommendations=recommendations
        )
