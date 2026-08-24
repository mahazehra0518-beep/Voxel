"""
VOXEL Base Models and Data Contracts.
Defines unified result structures, predictor statuses, and abstract base classes
for bioinformatics prediction adapters (VaxiJen, ToxinPred2, AllerTOP, etc.).
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any, List


class PredictorStatus(str, Enum):
    SUCCESS = "SUCCESS"
    FALLBACK_USED = "FALLBACK_USED"
    INPUT_REQUIRED = "INPUT_REQUIRED"
    ERROR = "ERROR"


class PredictionCategory(str, Enum):
    ANTIGENICITY = "Antigenicity"
    ALLERGENICITY = "Allergenicity"
    TOXICITY = "Toxicity"


@dataclass
class PredictionResult:
    """
    Standardized result contract for all biological prediction tools in VOXEL.
    Enforces scientific transparency, explicit thresholds, tool attribution,
    and computational disclaimers.
    """
    category: PredictionCategory
    tool_name: str
    prediction_score: Optional[float]
    classification: str
    threshold_used: float
    threshold_description: str
    interpretation: str
    status: PredictorStatus = PredictorStatus.SUCCESS
    details: Dict[str, Any] = field(default_factory=dict)
    disclaimer: str = (
        "Notice: Predictions are computational estimates derived from machine learning or "
        "bioinformatic models and require laboratory experimental validation."
    )

    def is_favorable_for_vaccine(self) -> bool:
        """
        Determines if the prediction meets the typical criteria for a vaccine candidate:
        - Antigenic: YES
        - Allergen: NO
        - Toxic: NO
        """
        cls_lower = self.classification.lower()
        if self.category == PredictionCategory.ANTIGENICITY:
            return "antigen" in cls_lower and "non" not in cls_lower
        elif self.category == PredictionCategory.ALLERGENICITY:
            return "non-allergen" in cls_lower or "non-allergenic" in cls_lower or "safe" in cls_lower
        elif self.category == PredictionCategory.TOXICITY:
            return "non-toxic" in cls_lower or "non-toxin" in cls_lower or "safe" in cls_lower
        return True


class BasePredictor:
    """
    Abstract base class for modular prediction adapters.
    """
    def __init__(self, name: str, category: PredictionCategory, default_threshold: float):
        self.name = name
        self.category = category
        self.default_threshold = default_threshold

    def predict(self, sequence: str, **kwargs) -> PredictionResult:
        raise NotImplementedError("Subclasses must implement the predict method.")
