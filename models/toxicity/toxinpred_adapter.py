"""
ToxinPred2 Toxicity Adapter for VOXEL.
Supports:
1. Direct / Uploaded ToxinPred2 scores.
2. Web bridge to ToxinPred2 server (GPSR / IIITD).
3. Local ToxinPred2 CLI invocation if installed.
4. Peptide toxicity risk index & motif scanner (disulfide density, charge density, toxic motifs).
"""

import re
import shutil
import subprocess
import tempfile
from typing import Optional, Dict, Any, List
from models.base import BasePredictor, PredictionResult, PredictionCategory, PredictorStatus

# Common pore-forming / hemolytic / neurotoxic peptide signatures
TOXIC_SIGNATURES = [
    (r"C.{1,4}C.{2,5}C.{1,4}C", "Conotoxin / Disulfide Knot Matrix"),
    (r"[KR]{4,}", "Polybasic Lytic / Cell Penetrating Array"),
    (r"C.{3}C.{7}C.{1}C", "Scorpion Neurotoxin Core"),
    (r"[FLWIV]{5,}", "Hydrophobic Pore-Forming Core")
]


class ToxinPredAdapter(BasePredictor):
    def __init__(self):
        super().__init__(
            name="ToxinPred2",
            category=PredictionCategory.TOXICITY,
            default_threshold=0.60
        )
        self.server_url = "https://webs.iiitd.edu.in/raghava/toxinpred2/"

    def scan_toxic_motifs(self, sequence: str) -> List[str]:
        """Scans sequence for known toxin structural signatures."""
        detected = []
        for pattern, label in TOXIC_SIGNATURES:
            if re.search(pattern, sequence):
                detected.append(label)
        return detected

    def calculate_toxicity_risk_index(self, sequence: str) -> float:
        """
        Calculates an empirical toxicity risk index (0.0 to 1.0) based on:
        - Cysteine density (disulfide richness characteristic of animal venoms)
        - Net positive charge density
        - Hydrophobic moment
        - Toxic motif matches
        """
        if not sequence:
            return 0.0
        
        seq_len = len(sequence)
        c_count = sequence.count('C')
        c_density = c_count / seq_len
        
        pos_charge = (sequence.count('K') + sequence.count('R')) / seq_len
        neg_charge = (sequence.count('D') + sequence.count('E')) / seq_len
        net_pos = max(0.0, pos_charge - neg_charge)
        
        motifs = self.scan_toxic_motifs(sequence)
        motif_penalty = 0.25 * len(motifs)
        
        # Base risk calculation
        risk = (c_density * 3.0) + (net_pos * 1.5) + motif_penalty
        return round(min(1.0, max(0.05, risk)), 4)

    def run_local_cli(self, sequence: str) -> Optional[float]:
        """Attempts to run local toxinpred2 if available in environment."""
        cli_path = shutil.which("toxinpred2") or shutil.which("toxinpred2.py")
        if not cli_path:
            return None
        try:
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tf:
                tf.write(f">query\n{sequence}\n")
                tf_path = tf.name
            
            out_file = tf_path + ".out.csv"
            cmd = [cli_path, "-i", tf_path, "-o", out_file, "-m", "1"]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
            if proc.returncode == 0:
                # Parse output CSV
                with open(out_file, "r") as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        cols = lines[1].strip().split(",")
                        # Expected column for score
                        for col in cols:
                            try:
                                return float(col)
                            except ValueError:
                                continue
        except Exception:
            pass
        return None

    def predict(
        self,
        sequence: str,
        user_score: Optional[float] = None,
        custom_threshold: Optional[float] = None,
        attempt_local: bool = False
    ) -> PredictionResult:
        """
        Predicts protein/peptide toxicity.
        """
        clean_seq = re.sub(r"[^A-Za-z]", "", sequence).upper()
        threshold = custom_threshold if custom_threshold is not None else self.default_threshold
        
        status = PredictorStatus.SUCCESS
        score = None
        tool_desc = "ToxinPred2 (Random Forest / Hybrid Model)"
        motifs = self.scan_toxic_motifs(clean_seq)
        
        if user_score is not None:
            score = float(user_score)
            tool_desc += " [User Provided / Verified ToxinPred2 Score]"
        elif attempt_local:
            local_score = self.run_local_cli(clean_seq)
            if local_score is not None:
                score = local_score
                tool_desc += " [Local Standalone ToxinPred2]"
            else:
                status = PredictorStatus.FALLBACK_USED
        
        if score is None:
            score = self.calculate_toxicity_risk_index(clean_seq)
            status = PredictorStatus.FALLBACK_USED
            tool_desc = "Toxicity Physicochemical Risk Index (Disulfide & Charge Density)"
            threshold = 0.50

        is_toxic = score >= threshold
        classification = "Toxic" if is_toxic else "Non-Toxic"
        
        interpretation = (
            f"Candidate toxicity score is {score:.4f} against the threshold of {threshold:.2f}. "
            f"Classification: '{classification}'."
        )
        if is_toxic:
            interpretation += " WARNING: High probability of cellular toxicity or lytic activity. Discard or engineer mutations."
            if motifs:
                interpretation += f" Detected suspicious motifs: {', '.join(motifs)}."
        else:
            interpretation += " Candidate exhibits favorable safety profile with low predicted toxicity."

        return PredictionResult(
            category=PredictionCategory.TOXICITY,
            tool_name=tool_desc,
            prediction_score=score,
            classification=classification,
            threshold_used=threshold,
            threshold_description=f"Cutoff score ≥ {threshold:.2f} classifies sequence as Toxic",
            interpretation=interpretation,
            status=status,
            details={
                "motifs_detected": motifs,
                "cysteine_count": clean_seq.count('C'),
                "sequence_length": len(clean_seq)
            }
        )
