"""Unit tests for src/tools/forcefield_dr2.py parsers and lookup logic."""

import pytest

from src.tools.forcefield_dr2 import (
    AtomicForcefield,
    ParamParser,
    PatchActionType,
    PatchParser,
    ResidueClassParser,
    ResidueClassRule,
    TopParser,
)


# ---------------------------------------------------------------------------
# ResidueClassRule.matches
# ---------------------------------------------------------------------------

def test_residue_class_rule_matches_all_aa():
    rule = ResidueClassRule("CTER", "all", present=["OXT"], absent=["H2"])
    assert rule.matches("ALA", ["N", "CA", "C", "O", "OXT"]) is True


def test_residue_class_rule_rejects_wrong_aa():
    rule = ResidueClassRule("CTER", ["GLY"], present=["OXT"], absent=[])
    assert rule.matches("ALA", ["OXT"]) is False


def test_residue_class_rule_rejects_absent_atom_present():
    rule = ResidueClassRule("CTER", "all", present=[], absent=["H2"])
    assert rule.matches("ALA", ["H2"]) is False


def test_residue_class_rule_rejects_missing_required_atom():
    rule = ResidueClassRule("CTER", "all", present=["OXT"], absent=[])
    assert rule.matches("ALA", ["N", "CA"]) is False


# ---------------------------------------------------------------------------
# PatchParser._parse_action_type
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "s,expected",
    [("MODIFY", PatchActionType.MODIFY), ("ADD", PatchActionType.ADD), ("DELETE", PatchActionType.DELETE)],
)
def test_parse_action_type(s, expected):
    assert PatchParser._parse_action_type(s) == expected


def test_parse_action_type_unknown_raises():
    with pytest.raises(ValueError):
        PatchParser._parse_action_type("BOGUS")


# ---------------------------------------------------------------------------
# ParamParser.parse
# ---------------------------------------------------------------------------

def test_param_parser_parses_nonbonded_line(tmp_path):
    param_file = tmp_path / "test.param"
    param_file.write_text("# comment\nNONBonded CT1 0.0200 2.2750 0.0100 1.9000\n\n")
    params = ParamParser.parse(str(param_file))
    assert set(params.keys()) == {"CT1"}
    p = params["CT1"]
    assert (p.epsilon_main, p.sigma_main, p.epsilon_14, p.sigma_14) == (0.02, 2.275, 0.01, 1.9)


def test_param_parser_unparsable_line_raises(tmp_path):
    param_file = tmp_path / "bad.param"
    param_file.write_text("GARBAGE LINE\n")
    with pytest.raises(ValueError):
        ParamParser.parse(str(param_file))


# ---------------------------------------------------------------------------
# TopParser.parse
# ---------------------------------------------------------------------------

def test_top_parser_parses_atom_line(tmp_path):
    top_file = tmp_path / "test.top"
    top_file.write_text("ALA atom CA TYPE=CT1 CHARGE=0.070 end\n")
    rows = TopParser.parse(str(top_file))
    assert len(rows) == 1
    row = rows[0]
    assert row.residue_name == "ALA"
    assert row.atom_name == "CA"
    assert row["type"] == "CT1"
    assert row["charge"] == 0.070


def test_top_parser_skips_comments_and_blank_lines(tmp_path):
    top_file = tmp_path / "test.top"
    top_file.write_text("# a comment\n\nALA atom CA TYPE=CT1 CHARGE=0.070 end\n")
    rows = TopParser.parse(str(top_file))
    assert len(rows) == 1


def test_top_parser_skips_unmatched_lines(tmp_path):
    top_file = tmp_path / "test.top"
    top_file.write_text("not a valid line at all\n")
    rows = TopParser.parse(str(top_file))
    assert rows == []


# ---------------------------------------------------------------------------
# ResidueClassParser.parse
# ---------------------------------------------------------------------------

def test_residue_class_parser_parses_header_and_atoms(tmp_path):
    rc_file = tmp_path / "residue-classes"
    rc_file.write_text("CTER: name=all present(OXT) absent(H2)\n")
    rules = ResidueClassParser.parse(str(rc_file))
    assert len(rules) == 1
    rule = rules[0]
    assert rule.class_name == "CTER"
    assert rule.amino_acid_names == "all"
    assert rule.present_atom_names == ["OXT"]
    assert rule.absent_atom_names == ["H2"]


def test_residue_class_parser_specific_aa(tmp_path):
    rc_file = tmp_path / "residue-classes"
    rc_file.write_text("NTER: name=GLY present(H1, H2, H3)\n")
    rules = ResidueClassParser.parse(str(rc_file))
    assert rules[0].amino_acid_names == ["GLY"]
    assert rules[0].present_atom_names == ["H1", "H2", "H3"]


# ---------------------------------------------------------------------------
# PatchParser.parse
# ---------------------------------------------------------------------------

def test_patch_parser_parses_string_and_numeric_vars(tmp_path):
    patch_file = tmp_path / "patch.top"
    patch_file.write_text("CTER MODIFY ATOM CA TYPE=CT2 CHARGE=-0.500\n")
    actions = PatchParser.parse(str(patch_file))
    assert len(actions) == 1
    action = actions[0]
    assert action.type == PatchActionType.MODIFY
    assert action.selection.residue_type == "CTER"
    assert action.selection.atom_name == "CA"
    assert action["TYPE"] == "CT2"
    assert action["CHARGE"] == -0.500


def test_patch_parser_skips_comments(tmp_path):
    patch_file = tmp_path / "patch.top"
    patch_file.write_text("# comment\n! another comment\nCTER MODIFY ATOM CA TYPE=CT2\n")
    actions = PatchParser.parse(str(patch_file))
    assert len(actions) == 1


def test_patch_parser_unmatched_line_raises(tmp_path):
    patch_file = tmp_path / "patch.top"
    patch_file.write_text("this is not a valid patch line\n")
    with pytest.raises(ValueError):
        PatchParser.parse(str(patch_file))


# ---------------------------------------------------------------------------
# AtomicForcefield (small 4-file fixture directory)
# ---------------------------------------------------------------------------

@pytest.fixture
def ff_dir(tmp_path):
    (tmp_path / "protein-allhdg5-5_new.top").write_text(
        "ALA atom CA TYPE=CT1 CHARGE=0.070 end\n"
    )
    (tmp_path / "protein-allhdg5-4_new.param").write_text(
        "NONBonded CT1 0.0200 2.2750 0.0100 1.9000\n"
        "NONBonded CT2 0.0300 2.0000 0.0200 1.8000\n"
    )
    (tmp_path / "residue-classes").write_text("CTER: name=all present(OXT) absent(H2)\n")
    (tmp_path / "patch.top").write_text("CTER MODIFY ATOM CA TYPE=CT2 CHARGE=-0.500\n")
    return str(tmp_path)


def test_atomic_forcefield_returns_base_vdw_without_patch(ff_dir):
    ff = AtomicForcefield(ff_dir)
    vdw = ff.get_vanderwaals("ALA", "CA", atom_list=[])  # no OXT -> no matching class
    assert (vdw.epsilon_main, vdw.sigma_main) == (0.02, 2.275)


def test_atomic_forcefield_applies_patch_override_for_vdw(ff_dir):
    ff = AtomicForcefield(ff_dir)
    vdw = ff.get_vanderwaals("ALA", "CA", atom_list=["OXT"])  # matches CTER -> patched to CT2
    assert (vdw.epsilon_main, vdw.sigma_main) == (0.03, 2.0)


def test_atomic_forcefield_returns_base_charge_without_patch(ff_dir):
    ff = AtomicForcefield(ff_dir)
    assert ff.get_charge("ALA", "CA", atom_list=[]) == 0.070


def test_atomic_forcefield_applies_patch_override_for_charge(ff_dir):
    ff = AtomicForcefield(ff_dir)
    assert ff.get_charge("ALA", "CA", atom_list=["OXT"]) == -0.500


def test_atomic_forcefield_unknown_atom_returns_zero_vdw(ff_dir):
    ff = AtomicForcefield(ff_dir)
    vdw = ff.get_vanderwaals("GLY", "N", atom_list=[])
    assert (vdw.epsilon_main, vdw.sigma_main, vdw.epsilon_14, vdw.sigma_14) == (0, 0, 0, 0)
