"""
Module 2: Epitope Analysis, Multi-Criteria Filtering, Exclusion, and Prioritization Engine for VOXEL.
Supports:
1. Flexible dataset ingestion (CSV, TSV, XLSX) and automated column mapping.
2. Distinct 3-tier criteria architecture:
   - Filtering Criteria (IC50 affinity, binding rank, antigenicity).
   - Exclusion Criteria (Toxicity, Allergenicity, Ambiguous AAs).
   - Prioritization Criteria (Normalized MCDA Composite Scoring & Promiscuity).
3. Comprehensive Rejection Reason Logger and Auditing Trail.
4. Summary analytics and prioritized list generation.
"""

import math
import re
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any


@dataclass
class EpitopeFilteringConfig:
    # Tier 1: Filtering Criteria (Range / Thresholds)
    enable_ic50_filter: bool = True
    max_ic50_nm: float = 500.0  # Common benchmark: <50 nM (Strong), <500 nM (Intermediate)
    
    enable_rank_filter: bool = True
    max_binding_rank: float = 2.0  # Percentile rank (e.g. <= 0.5% strong, <= 2.0% weak)
    
    enable_antigenicity_filter: bool = False
    min_antigenicity_score: float = 0.50
    
    # Tier 2: Exclusion Criteria (Safety / Quality Filters)
    exclude_toxic: bool = True
    toxin_score_cutoff: float = 0.60
    
    exclude_allergenic: bool = True
    allergen_score_cutoff: float = 0.50
    
    exclude_ambiguous_aa: bool = True
    
    # Tier 3: Prioritization Weights (MCDA - Multi-Criteria Decision Analysis)
    weight_ic50: float = 0.35
    weight_rank: float = 0.25
    weight_antigenicity: float = 0.25
    weight_promiscuity: float = 0.15


@dataclass
class CandidateEvaluation:
    index: int
    epitope: str
    allele: str
    passed: bool
    rejection_reasons: List[str] = field(default_factory=list)
    priority_score: Optional[float] = None
    binding_category: str = "Unknown"


@dataclass
class EpitopeAnalysisResult:
    total_candidates: int
    passed_candidates_count: int
    rejected_candidates_count: int
    rejection_breakdown: Dict[str, int]
    prioritized_df: pd.DataFrame
    audited_full_df: pd.DataFrame
    column_mapping: Dict[str, str]
    config: EpitopeFilteringConfig


class EpitopeDatasetEngine:
    # Common column name synonyms in immunoinformatics tools (IEDB, NetMHCpan, Syfpeithi, etc.)
    SYNONYMS = {
        "epitope": ["epitope", "peptide", "sequence", "core", "pep_seq", "peptide_sequence", "pep", "mer"],
        "allele": ["allele", "mhc", "hla", "mhc_allele", "gene", "hla_allele", "molecule"],
        "ic50": ["ic50", "ic50_nm", "affinity", "binding_affinity", "ic50(nm)", "affinity(nm)", "score_nm"],
        "rank": ["rank", "%rank", "percentile_rank", "percentile", "binding_rank", "score_rank", "el_rank"],
        "antigenicity": ["antigenicity", "vaxijen", "vaxijen_score", "antigenic_score", "protective_score", "ag_score"],
        "toxicity": ["toxicity", "toxin", "toxinpred", "is_toxic", "toxic", "toxin_score", "tx_score"],
        "allergenicity": ["allergenicity", "allertop", "allergen", "is_allergen", "allergen_score", "al_score"]
    }

    @classmethod
    def auto_detect_columns(cls, df: pd.DataFrame) -> Dict[str, Optional[str]]:
        """
        Scans DataFrame columns and maps them to standard fields.
        """
        clean_cols = {col: re.sub(r"[^a-zA-Z0-9%]", "", str(col).lower()) for col in df.columns}
        mapped: Dict[str, Optional[str]] = {
            "epitope": None,
            "allele": None,
            "ic50": None,
            "rank": None,
            "antigenicity": None,
            "toxicity": None,
            "allergenicity": None
        }

        for std_key, syn_list in cls.SYNONYMS.items():
            for col, norm_col in clean_cols.items():
                if mapped[std_key] is not None:
                    break
                for syn in syn_list:
                    norm_syn = re.sub(r"[^a-zA-Z0-9%]", "", syn.lower())
                    if norm_col == norm_syn or norm_syn in norm_col:
                        mapped[std_key] = col
                        break

        return mapped

    @classmethod
    def parse_dataset(cls, file_obj, filename: str) -> pd.DataFrame:
        """
        Parses uploaded CSV, TSV, or Excel file into a pandas DataFrame.
        """
        fn_lower = filename.lower()
        if fn_lower.endswith(".csv"):
            # Auto-detect comma or semicolon
            df = pd.read_csv(file_obj, sep=None, engine="python")
        elif fn_lower.endswith(".tsv") or fn_lower.endswith(".txt"):
            df = pd.read_csv(file_obj, sep="\t")
        elif fn_lower.endswith(".xlsx") or fn_lower.endswith(".xls"):
            df = pd.read_excel(file_obj)
        else:
            # Fallback to general read_csv
            df = pd.read_csv(file_obj)
        
        # Clean column names (strip leading/trailing whitespaces)
        df.columns = [str(c).strip() for c in df.columns]
        return df

    @classmethod
    def is_value_toxic(cls, val: Any, cutoff: float = 0.60) -> Tuple[bool, Optional[float]]:
        """Determines if a value represents a toxic candidate."""
        if pd.isna(val) or val is None or str(val).strip() == "":
            return False, None
        
        val_str = str(val).strip().lower()
        if val_str in ["toxic", "yes", "true", "positive", "1", "t"]:
            return True, 1.0
        if val_str in ["non-toxic", "nontoxic", "no", "false", "negative", "0", "nt", "safe"]:
            return False, 0.0
        
        try:
            num = float(val)
            return num >= cutoff, num
        except ValueError:
            return "toxic" in val_str and "non" not in val_str, None

    @classmethod
    def is_value_allergenic(cls, val: Any, cutoff: float = 0.50) -> Tuple[bool, Optional[float]]:
        """Determines if a value represents an allergenic candidate."""
        if pd.isna(val) or val is None or str(val).strip() == "":
            return False, None
        
        val_str = str(val).strip().lower()
        if val_str in ["allergen", "allergenic", "yes", "true", "positive", "1", "al"]:
            return True, 1.0
        if val_str in ["non-allergen", "non-allergenic", "nonallergen", "no", "false", "negative", "0", "safe"]:
            return False, 0.0
            
        try:
            num = float(val)
            return num >= cutoff, num
        except ValueError:
            return "allergen" in val_str and "non" not in val_str, None

    @classmethod
    def evaluate_dataset(
        cls,
        df: pd.DataFrame,
        mapping: Dict[str, Optional[str]],
        config: EpitopeFilteringConfig
    ) -> EpitopeAnalysisResult:
        """
        Executes complete multi-criteria screening, exclusion, rejection logging,
        and MCDA prioritization across all dataset rows.
        """
        audited_df = df.copy()
        
        epitope_col = mapping.get("epitope")
        allele_col = mapping.get("allele")
        ic50_col = mapping.get("ic50")
        rank_col = mapping.get("rank")
        ag_col = mapping.get("antigenicity")
        tx_col = mapping.get("toxicity")
        al_col = mapping.get("allergenicity")

        total_rows = len(audited_df)
        rejection_reasons_list: List[List[str]] = []
        pass_flags: List[bool] = []
        priority_scores: List[Optional[float]] = []
        binding_categories: List[str] = []
        
        rejection_breakdown = {
            "IC50 Binding Threshold": 0,
            "Binding Percentile Rank": 0,
            "Low Antigenicity": 0,
            "Excluded: Toxic": 0,
            "Excluded: Allergenic": 0,
            "Non-Standard / Empty Sequence": 0
        }

        # Calculate allele promiscuity per unique epitope across the dataset
        allele_counts_per_pep: Dict[str, int] = {}
        if epitope_col and epitope_col in audited_df.columns:
            if allele_col and allele_col in audited_df.columns:
                # Count distinct alleles for each clean peptide
                grouped = audited_df.groupby(audited_df[epitope_col].astype(str).str.strip().str.upper())[allele_col].nunique()
                allele_counts_per_pep = grouped.to_dict()
            else:
                allele_counts_per_pep = {str(k).strip().upper(): 1 for k in audited_df[epitope_col].unique()}

        max_promiscuity = max(allele_counts_per_pep.values()) if allele_counts_per_pep else 1

        for idx, row in audited_df.iterrows():
            row_reasons: List[str] = []
            
            # --- 1. Sequence Validation ---
            pep_raw = str(row[epitope_col]).strip().upper() if epitope_col and pd.notna(row[epitope_col]) else ""
            pep_clean = re.sub(r"[^A-Z]", "", pep_raw)
            
            if not pep_clean or len(pep_clean) < 7:
                row_reasons.append(f"Invalid sequence length ({len(pep_clean)} AA). Must be ≥ 7 AA.")
                rejection_breakdown["Non-Standard / Empty Sequence"] += 1
            elif config.exclude_ambiguous_aa and set(pep_clean) - set("ACDEFGHIKLMNPQRSTVWY"):
                diff = set(pep_clean) - set("ACDEFGHIKLMNPQRSTVWY")
                row_reasons.append(f"Contains non-standard amino acids: {', '.join(sorted(diff))}")
                rejection_breakdown["Non-Standard / Empty Sequence"] += 1

            # --- 2. IC50 Binding Filtering Criteria ---
            ic50_val: Optional[float] = None
            if ic50_col and pd.notna(row[ic50_col]):
                try:
                    ic50_val = float(row[ic50_col])
                    if config.enable_ic50_filter and ic50_val > config.max_ic50_nm:
                        row_reasons.append(f"IC50 affinity ({ic50_val:.1f} nM) exceeds cutoff ({config.max_ic50_nm:.1f} nM)")
                        rejection_breakdown["IC50 Binding Threshold"] += 1
                except ValueError:
                    pass

            # Binding Category classification
            if ic50_val is not None:
                if ic50_val < 50.0:
                    b_cat = "Strong Binder (<50 nM)"
                elif ic50_val < 500.0:
                    b_cat = "Intermediate Binder (<500 nM)"
                elif ic50_val < 5000.0:
                    b_cat = "Weak Binder (<5000 nM)"
                else:
                    b_cat = "Non-Binder (≥5000 nM)"
            else:
                b_cat = "Not Specified"
            binding_categories.append(b_cat)

            # --- 3. Percentile Rank Filtering Criteria ---
            rank_val: Optional[float] = None
            if rank_col and pd.notna(row[rank_col]):
                try:
                    rank_val = float(str(row[rank_col]).replace("%", "").strip())
                    if config.enable_rank_filter and rank_val > config.max_binding_rank:
                        row_reasons.append(f"Binding rank ({rank_val:.2f}%) exceeds threshold ({config.max_binding_rank:.2f}%)")
                        rejection_breakdown["Binding Percentile Rank"] += 1
                except ValueError:
                    pass

            # --- 4. Antigenicity Filtering Criteria ---
            ag_val: Optional[float] = None
            if ag_col and pd.notna(row[ag_col]):
                try:
                    ag_val = float(row[ag_col])
                    if config.enable_antigenicity_filter and ag_val < config.min_antigenicity_score:
                        row_reasons.append(f"Antigenicity score ({ag_val:.3f}) below threshold ({config.min_antigenicity_score:.2f})")
                        rejection_breakdown["Low Antigenicity"] += 1
                except ValueError:
                    pass

            # --- 5. Toxicity Exclusion Criteria ---
            if tx_col and config.exclude_toxic and pd.notna(row[tx_col]):
                is_tx, tx_num = cls.is_value_toxic(row[tx_col], cutoff=config.toxin_score_cutoff)
                if is_tx:
                    score_info = f" (Score: {tx_num:.2f})" if tx_num is not None else ""
                    row_reasons.append(f"Excluded: Flagged as Toxic{score_info}")
                    rejection_breakdown["Excluded: Toxic"] += 1

            # --- 6. Allergenicity Exclusion Criteria ---
            if al_col and config.exclude_allergenic and pd.notna(row[al_col]):
                is_al, al_num = cls.is_value_allergenic(row[al_col], cutoff=config.allergen_score_cutoff)
                if is_al:
                    score_info = f" (Score: {al_num:.2f})" if al_num is not None else ""
                    row_reasons.append(f"Excluded: Flagged as Allergenic{score_info}")
                    rejection_breakdown["Excluded: Allergenic"] += 1

            passed = (len(row_reasons) == 0)
            pass_flags.append(passed)
            rejection_reasons_list.append(row_reasons)

            # --- 7. Prioritization MCDA Scoring (for passing and full ranking) ---
            # Normalized components:
            # Affinity score (higher is better, 0 to 1) using logarithmic affinity scale
            if ic50_val is not None and ic50_val > 0:
                s_affinity = max(0.0, min(1.0, 1.0 - (math.log10(min(50000.0, max(1.0, ic50_val))) / math.log10(50000.0))))
            else:
                s_affinity = 0.5

            # Rank score (lower %rank is better, 0 to 1)
            if rank_val is not None:
                s_rank = max(0.0, min(1.0, 1.0 - (min(10.0, rank_val) / 10.0)))
            else:
                s_rank = 0.5

            # Antigenicity score (0 to 1)
            if ag_val is not None:
                s_ag = max(0.0, min(1.0, ag_val))
            else:
                s_ag = 0.5

            # Promiscuity score (0 to 1)
            pep_promiscuity = allele_counts_per_pep.get(pep_clean, 1)
            s_prom = min(1.0, pep_promiscuity / max(1, max_promiscuity))

            # Total weighted priority score (0 to 100)
            w_total = config.weight_ic50 + config.weight_rank + config.weight_antigenicity + config.weight_promiscuity
            if w_total > 0:
                composite = (
                    (config.weight_ic50 * s_affinity) +
                    (config.weight_rank * s_rank) +
                    (config.weight_antigenicity * s_ag) +
                    (config.weight_promiscuity * s_prom)
                ) / w_total
                priority_scores.append(round(composite * 100.0, 2))
            else:
                priority_scores.append(50.0)

        # Append audit columns to DataFrame
        audited_df["Screening_Status"] = ["PASSED" if p else "REJECTED" for p in pass_flags]
        audited_df["Rejection_Reason"] = [
            "PASSED (Meets all screening criteria)" if len(r) == 0 else " ; ".join(r)
            for r in rejection_reasons_list
        ]
        audited_df["Binding_Affinity_Category"] = binding_categories
        audited_df["Priority_Score"] = priority_scores
        
        # Add Promiscuity (Distinct Allele Count)
        if epitope_col:
            audited_df["Allele_Promiscuity_Count"] = [
                allele_counts_per_pep.get(str(p).strip().upper(), 1) for p in audited_df[epitope_col]
            ]

        # Prioritized DataFrame (Filtered to PASSED candidates only, sorted by priority score)
        prioritized_df = audited_df[audited_df["Screening_Status"] == "PASSED"].copy()
        
        # Sort by Priority_Score descending, then by IC50 ascending if available
        sort_cols = ["Priority_Score"]
        ascending_flags = [False]
        if ic50_col and ic50_col in prioritized_df.columns:
            sort_cols.append(ic50_col)
            ascending_flags.append(True)
        if rank_col and rank_col in prioritized_df.columns:
            sort_cols.append(rank_col)
            ascending_flags.append(True)

        prioritized_df = prioritized_df.sort_values(by=sort_cols, ascending=ascending_flags).reset_index(drop=True)
        prioritized_df.insert(0, "Priority_Rank", range(1, len(prioritized_df) + 1))

        passed_count = sum(pass_flags)
        rejected_count = total_rows - passed_count

        return EpitopeAnalysisResult(
            total_candidates=total_rows,
            passed_candidates_count=passed_count,
            rejected_candidates_count=rejected_count,
            rejection_breakdown=rejection_breakdown,
            prioritized_df=prioritized_df,
            audited_full_df=audited_df,
            column_mapping=mapping,
            config=config
        )
