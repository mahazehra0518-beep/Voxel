"""
VaxiJen Antigenicity Adapter for VOXEL.
Supports:
1. Target organism-specific thresholds (Bacteria: 0.4, Virus: 0.4, Tumor: 0.5, Parasite: 0.5, Fungal: 0.5).
2. Remote server bridge with graceful timeout handling.
3. User-supplied VaxiJen prediction scores and uploaded results.
4. Physicochemical Antigenic Propensity Index based on the Kolaskar & Tongaonkar (1990) scale.
"""

import re
import requests
from typing import Optional, Dict, Any
from models.base import BasePredictor, PredictionResult, PredictionCategory, PredictorStatus

# Kolaskar & Tongaonkar (1990) semi-empirical antigenic propensity scale
KOLASKAR_TONGAONKAR_SCALE = {
    'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
    'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
    'M': 0.826, 'N': 0.776, 'P': 1.064, 'Q': 1.015, 'R': 0.971,
    'S': 1.009, 'T': 0.909, 'V': 1.383, 'W': 0.893, 'Y': 1.161
}

DEFAULT_ORGANISM_THRESHOLDS = {
    "Bacteria": 0.4,
    "Virus": 0.4,
    "Tumour": 0.5,
    "Parasite": 0.5,
    "Fungi": 0.5,
    "General/Default": 0.4
}


class VaxiJenAdapter(BasePredictor):
    def __init__(self):
        super().__init__(
            name="VaxiJen v2.0",
            category=PredictionCategory.ANTIGENICITY,
            default_threshold=0.4
        )
        self.server_url = "http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html"

    def calculate_kolaskar_propensity(self, sequence: str) -> float:
        """Calculates mean Kolaskar-Tongaonkar antigenic propensity for the sequence."""
        clean_seq = "".join([c.upper() for c in sequence if c.upper() in KOLASKAR_TONGAONKAR_SCALE])
        if not clean_seq:
            return 1.0
        scores = [KOLASKAR_TONGAONKAR_SCALE[aa] for aa in clean_seq]
        return round(sum(scores) / len(scores), 4)

    def query_remote_server(self, sequence: str, organism: str = "Virus", timeout: float = 8.0) -> Optional[float]:
        """Attempts to query the public VaxiJen v2.0 web server."""
        try:
            # Map organism to form parameter values expected by VaxiJen
            target_map = {
                "Bacteria": "bacteria",
                "Virus": "virus",
                "Tumour": "tumor",
                "Parasite": "parasite",
                "Fungi": "fungal"
            }
            target_val = target_map.get(organism, "virus")
            
            # Post FASTA payload
            payload = {
                "seq": f">seq\n{sequence}",
                "Target": target_val
            }
            headers = {"User-Agent": "VOXEL-Platform/1.0 (Bioinformatics Research Client)"}
            resp = requests.post(self.server_url, data=payload, headers=headers, timeout=timeout)
            
            if resp.status_code == 200:
                # Search for overall prediction score in HTML response
                # Typical pattern: "Overall Prediction Protective Antigen = 0.5432"
                match = re.search(r"Overall Prediction[^=]*=\s*([0-9]+\.[0-9]+)", resp.text, re.IGNORECASE)
                if match:
                    return float(match.group(1))
        except Exception:
            pass
        return None

    def predict(
        self,
        sequence: str,
        user_score: Optional[float] = None,
        organism: str = "Virus",
        custom_threshold: Optional[float] = None,
        attempt_remote: bool = False
    ) -> PredictionResult:
        """
        Evaluates antigenicity using provided score, remote server bridge, or fallback ACC/propensity index.
        """
        clean_seq = re.sub(r"[^A-Za-z]", "", sequence).upper()
        threshold = custom_threshold if custom_threshold is not None else DEFAULT_ORGANISM_THRESHOLDS.get(organism, 0.4)
        
        status = PredictorStatus.SUCCESS
        score = None
        tool_desc = f"VaxiJen v2.0 (Organism model: {organism})"
        
        if user_score is not None:
            score = float(user_score)
            status = PredictorStatus.SUCCESS
            tool_desc += " [User Provided / Verified VaxiJen Score]"
        elif attempt_remote:
            remote_score = self.query_remote_server(clean_seq, organism=organism)
            if remote_score is not None:
                score = remote_score
                status = PredictorStatus.SUCCESS
                tool_desc += " [Live Web Server Bridge]"
            else:
                status = PredictorStatus.FALLBACK_USED
        
        # If no score was obtained from user or remote server, compute Kolaskar-Tongaonkar propensity
        kt_score = self.calculate_kolaskar_propensity(clean_seq)
        
        if score is None:
            # Normalized approximation based on Kolaskar-Tongaonkar mean propensity (baseline ~1.0)
            # A score > 1.0 indicates higher than average antigenic propensity
            score = round((kt_score - 0.7) / 0.6, 4)
            score = max(0.0, min(1.0, score))
            status = PredictorStatus.FALLBACK_USED
            tool_desc = f"Kolaskar-Tongaonkar Antigenic Index (Scale Mean: {kt_score})"
            threshold = 0.5

        is_antigen = score >= threshold
        classification = "Probable Antigen" if is_antigen else "Probable Non-Antigen"
        
        interpretation = (
            f"Candidate has an antigenicity score of {score:.4f} against the configured threshold of "
            f"{threshold:.2f} ({organism} model). It is classified as '{classification}'."
        )
        if not is_antigen:
            interpretation += " Low antigenicity indicates reduced likelihood of inducing protective immune response."
        else:
            interpretation += " Favorable for triggering humoral/cellular immune recognition."

        return PredictionResult(
            category=PredictionCategory.ANTIGENICITY,
            tool_name=tool_desc,
            prediction_score=score,
            classification=classification,
            threshold_used=threshold,
            threshold_description=f"Cutoff score ≥ {threshold:.2f} for {organism} indicates protective antigen",
            interpretation=interpretation,
            status=status,
            details={
                "organism": organism,
                "kolaskar_tongaonkar_mean": kt_score,
                "sequence_length": len(clean_seq)
            }
        )
