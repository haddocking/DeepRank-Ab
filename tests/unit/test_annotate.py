"""Unit tests for pure/near-pure helpers in src/tools/annotate.py."""

import pytest

from src.tools.annotate import (
    _norm_role_from_id,
    extract_chain_sequence,
    locate,
    region_h,
    region_l,
)


# ---------------------------------------------------------------------------
# region_l / region_h
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "position,expected",
    [
        (26, "CONST"),
        (27, "L1"),
        (38, "L1"),
        (39, "CONST"),
        (55, "CONST"),
        (56, "L2"),
        (65, "L2"),
        (66, "CONST"),
        (104, "CONST"),
        (105, "L3"),
        (117, "L3"),
        (118, "CONST"),
    ],
)
def test_region_l_boundaries(position, expected):
    assert region_l(position) == expected


@pytest.mark.parametrize(
    "position,expected",
    [
        (26, "CONST"),
        (27, "H1"),
        (38, "H1"),
        (39, "CONST"),
        (56, "H2"),
        (65, "H2"),
        (105, "H3"),
        (117, "H3"),
        (118, "CONST"),
    ],
)
def test_region_h_boundaries(position, expected):
    assert region_h(position) == expected


# ---------------------------------------------------------------------------
# _norm_role_from_id
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "rec_id,expected",
    [
        ("foo.H", "H"),
        ("foo.L", "L"),
        ("model_001.h", "H"),
        ("  foo.L  ", "L"),
        ("", ""),
    ],
)
def test_norm_role_from_id(rec_id, expected):
    assert _norm_role_from_id(rec_id) == expected


# ---------------------------------------------------------------------------
# extract_chain_sequence
# ---------------------------------------------------------------------------

class _StubResidue:
    def __init__(self, resname):
        self._resname = resname

    def get_resname(self):
        return self._resname


def test_extract_chain_sequence_maps_known_residues():
    chain = [_StubResidue("ALA"), _StubResidue("GLY"), _StubResidue("MSE")]
    assert extract_chain_sequence(chain) == "AGM"


def test_extract_chain_sequence_unknown_residue_is_x():
    chain = [_StubResidue("ALA"), _StubResidue("HOH")]
    assert extract_chain_sequence(chain) == "AX"


def test_extract_chain_sequence_empty_chain():
    assert extract_chain_sequence([]) == ""


# ---------------------------------------------------------------------------
# locate
# ---------------------------------------------------------------------------

def test_locate_finds_exact_substring():
    fullseq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV"
    subseq = "QISFVKSH"
    start, end = locate(subseq, fullseq, tag="test")
    assert fullseq[start:end] == subseq


def test_locate_raises_on_no_alignment():
    with pytest.raises(RuntimeError):
        locate("", "", tag="test")
