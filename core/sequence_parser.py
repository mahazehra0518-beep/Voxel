"""
Sequence Parser and Physicochemical Analyzer for VOXEL.
Uses Biopython (Bio.SeqUtils.ProtParam) for rigorous protein sequence parsing,
validation, and molecular characterization.
"""

import io
import re
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Standard IUPAC 20 amino acids
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class PhysicochemicalProperties:
    length: int
    molecular_weight_da: float
    molecular_weight_kda: float
    isoelectric_point: float
    gravy_hydropathy: float
    instability_index: float
    is_stable_in_vitro: bool
    aromaticity: float
    secondary_structure_fractions: Dict[str, float]  # Helix, Turn, Sheet
    amino_acid_counts: Dict[str, int]
    amino_acid_percentages: Dict[str, float]
    charge_at_ph7: float


@dataclass
class ProteinSequence:
    identifier: str
    description: str
    raw_sequence: str
    clean_sequence: str
    organism: Optional[str] = None
    is_valid: bool = True
    validation_warnings: List[str] = field(default_factory=list)
    properties: Optional[PhysicochemicalProperties] = None


class SequenceParser:
    @staticmethod
    def clean_raw_sequence(seq_str: str) -> str:
        """Removes whitespaces, line breaks, numbers, and hyphens."""
        return re.sub(r"[\s\d\-_*]", "", seq_str).upper()

    @classmethod
    def parse_fasta(cls, fasta_content: str, default_organism: Optional[str] = None) -> List[ProteinSequence]:
        """
        Parses FASTA formatted string containing one or multiple records.
        If plain sequence text without headers is passed, handles it gracefully.
        """
        records: List[ProteinSequence] = []
        trimmed = fasta_content.strip()
        if not trimmed:
            return []

        # Check if text starts with FASTA header '>'
        if not trimmed.startswith(">"):
            # Plain sequence input
            clean_seq = cls.clean_raw_sequence(trimmed)
            record = cls._build_protein_sequence(
                seq_id="Candidate_Protein",
                description="User input sequence",
                raw_seq=clean_seq,
                organism=default_organism
            )
            return [record]

        # Use Bio.SeqIO to parse FASTA
        fasta_io = io.StringIO(trimmed)
        for record in SeqIO.parse(fasta_io, "fasta"):
            clean_seq = cls.clean_raw_sequence(str(record.seq))
            
            # Extract organism from description if present (e.g. [Homo sapiens] or OS=...)
            org = default_organism
            desc = record.description
            os_match = re.search(r"OS=([^=]+?)(?:\s+[A-Z]{2}=|$)", desc)
            bracket_match = re.search(r"\[([^\]]+)\]", desc)
            if os_match:
                org = os_match.group(1).strip()
            elif bracket_match:
                org = bracket_match.group(1).strip()

            prot_record = cls._build_protein_sequence(
                seq_id=record.id,
                description=desc,
                raw_seq=clean_seq,
                organism=org
            )
            records.append(prot_record)

        return records

    @classmethod
    def _build_protein_sequence(
        cls,
        seq_id: str,
        description: str,
        raw_seq: str,
        organism: Optional[str] = None
    ) -> ProteinSequence:
        warnings: List[str] = []
        is_valid = True

        if not raw_seq:
            return ProteinSequence(
                identifier=seq_id,
                description=description,
                raw_sequence="",
                clean_sequence="",
                organism=organism,
                is_valid=False,
                validation_warnings=["Sequence is empty."]
            )

        # Check for non-standard amino acid characters
        invalid_chars = set(raw_seq) - STANDARD_AMINO_ACIDS
        if invalid_chars:
            warnings.append(
                f"Contains non-standard/ambiguous amino acid codes: {', '.join(sorted(invalid_chars))}. "
                f"These were filtered for physicochemical parameter analysis."
            )
            # Filter non-standard characters for Biopython ProtParam
            seq_for_analysis = "".join([c for c in raw_seq if c in STANDARD_AMINO_ACIDS])
        else:
            seq_for_analysis = raw_seq

        if len(seq_for_analysis) < 5:
            warnings.append("Sequence length is under 5 amino acids. Some parameters cannot be reliably computed.")
            is_valid = False

        properties = None
        if len(seq_for_analysis) >= 5:
            try:
                analysis = ProteinAnalysis(seq_for_analysis)
                mw = analysis.molecular_weight()
                pi = analysis.isoelectric_point()
                gravy = analysis.gravy()
                instability = analysis.instability_index()
                aromaticity = analysis.aromaticity()
                sec_struct = analysis.secondary_structure_fraction()  # (Helix, Turn, Sheet)
                aa_counts = analysis.count_amino_acids()
                total_len = len(seq_for_analysis)
                aa_percentages = {aa: round((cnt / total_len) * 100, 2) for aa, cnt in aa_counts.items()}
                charge_ph7 = analysis.charge_at_pH(7.0)

                properties = PhysicochemicalProperties(
                    length=len(raw_seq),
                    molecular_weight_da=round(mw, 2),
                    molecular_weight_kda=round(mw / 1000.0, 3),
                    isoelectric_point=round(pi, 2),
                    gravy_hydropathy=round(gravy, 3),
                    instability_index=round(instability, 2),
                    is_stable_in_vitro=instability < 40.0,
                    aromaticity=round(aromaticity, 3),
                    secondary_structure_fractions={
                        "Helix": round(sec_struct[0] * 100, 1),
                        "Turn": round(sec_struct[1] * 100, 1),
                        "Sheet": round(sec_struct[2] * 100, 1),
                    },
                    amino_acid_counts=aa_counts,
                    amino_acid_percentages=aa_percentages,
                    charge_at_ph7=round(charge_ph7, 2)
                )
            except Exception as e:
                warnings.append(f"ProtParam analysis notice: {str(e)}")

        return ProteinSequence(
            identifier=seq_id,
            description=description,
            raw_sequence=raw_seq,
            clean_sequence=raw_seq,
            organism=organism,
            is_valid=is_valid,
            validation_warnings=warnings,
            properties=properties
        )
