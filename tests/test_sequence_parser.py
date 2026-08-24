"""
Unit tests for SequenceParser and ProteinSequence in VOXEL.
"""

import pytest
from core.sequence_parser import SequenceParser, ProteinSequence


def test_parse_fasta_header():
    fasta_text = ">test_protein [Homo sapiens]\nACDEFGHIKLMNPQRSTVWY\n"
    records = SequenceParser.parse_fasta(fasta_text)
    assert len(records) == 1
    rec = records[0]
    assert rec.identifier == "test_protein"
    assert rec.organism == "Homo sapiens"
    assert rec.clean_sequence == "ACDEFGHIKLMNPQRSTVWY"
    assert rec.is_valid is True
    assert rec.properties is not None
    assert rec.properties.length == 20
    assert rec.properties.molecular_weight_da > 0
    assert 0 < rec.properties.isoelectric_point < 14


def test_parse_plain_sequence():
    plain_seq = "  mrgshhhhhh gmasmtggqq mgrdlydddd kdrwgs  \n"
    records = SequenceParser.parse_fasta(plain_seq, default_organism="Synthetic")
    assert len(records) == 1
    rec = records[0]
    assert rec.identifier == "Candidate_Protein"
    assert rec.clean_sequence == "MRGSHHHHHHGMASMTGGQQMGRDLYDDDDKDRWGS"
    assert rec.properties is not None
    assert rec.properties.length == 36


def test_ambiguous_amino_acids_warning():
    seq_with_x = ">prot_x\nACDEFGHIKLMNPQRSTVWYXBZ\n"
    records = SequenceParser.parse_fasta(seq_with_x)
    assert len(records) == 1
    rec = records[0]
    assert rec.is_valid is True
    assert len(rec.validation_warnings) > 0
    assert "non-standard" in rec.validation_warnings[0]
    assert rec.properties is not None
    assert rec.properties.length == 23


def test_empty_sequence():
    records = SequenceParser.parse_fasta("")
    assert len(records) == 0

    records_empty_header = SequenceParser.parse_fasta(">empty_header\n\n")
    assert len(records_empty_header) == 1
    assert records_empty_header[0].is_valid is False
