"""Unit tests for pure/near-pure helpers in scripts/inference.py."""

import json

import h5py
import numpy as np
import pandas as pd
import pytest
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from scripts.inference import (
    _free_chain_id,
    _norm_chain,
    _read_fasta_sequences,
    ca_coords,
    clamp_cores,
    correct_json,
    flag_low_contacts,
    get_chain_sequence,
    hdf5_to_csv,
    three_to_one,
)


# ---------------------------------------------------------------------------
# _norm_chain
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw", ["-", "none", "None", "null", "NULL", "na", "", "  ", None])
def test_norm_chain_missing_values(raw):
    assert _norm_chain(raw) is None


@pytest.mark.parametrize("raw,expected", [("H", "H"), (" H ", "H"), ("L", "L")])
def test_norm_chain_passthrough(raw, expected):
    assert _norm_chain(raw) == expected


# ---------------------------------------------------------------------------
# _free_chain_id
# ---------------------------------------------------------------------------

def test_free_chain_id_skips_occupied_and_b():
    occupied = {"A", "B", "C"}
    assert _free_chain_id(occupied) == "D"


def test_free_chain_id_never_returns_b():
    occupied = set()
    assert _free_chain_id(occupied) == "A"
    occupied.add("A")
    assert _free_chain_id(occupied) == "C"


def test_free_chain_id_raises_when_exhausted():
    occupied = set("ACDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    with pytest.raises(RuntimeError):
        _free_chain_id(occupied)


# ---------------------------------------------------------------------------
# clamp_cores
# ---------------------------------------------------------------------------

def test_clamp_cores_respects_cpu_count(monkeypatch):
    monkeypatch.setattr("scripts.inference.os.cpu_count", lambda: 4)
    assert clamp_cores(requested=100, task_count=100) == 4


def test_clamp_cores_respects_task_count(monkeypatch):
    monkeypatch.setattr("scripts.inference.os.cpu_count", lambda: 100)
    assert clamp_cores(requested=100, task_count=3) == 3


def test_clamp_cores_never_below_one(monkeypatch):
    monkeypatch.setattr("scripts.inference.os.cpu_count", lambda: 8)
    assert clamp_cores(requested=0, task_count=0) == 1


# ---------------------------------------------------------------------------
# three_to_one / get_chain_sequence
# ---------------------------------------------------------------------------

def test_three_to_one_known_residues():
    mapping = three_to_one()
    assert mapping["ALA"] == "A"
    assert mapping["GLY"] == "G"
    assert mapping["TRP"] == "W"


def _build_structure(chain_id, resnames):
    structure = Structure("s")
    model = Model(0)
    structure.add(model)
    chain = Chain(chain_id)
    model.add(chain)
    for i, resname in enumerate(resnames, start=1):
        residue = Residue((" ", i, " "), resname, "")
        chain.add(residue)
    return structure


def test_get_chain_sequence_maps_known_residues():
    structure = _build_structure("H", ["ALA", "GLY", "TRP"])
    assert get_chain_sequence(structure, "H") == "AGW"


def test_get_chain_sequence_unknown_residue_is_x():
    structure = _build_structure("H", ["ALA", "XYZ"])
    assert get_chain_sequence(structure, "H") == "AX"


def test_get_chain_sequence_missing_chain_returns_empty():
    structure = _build_structure("H", ["ALA"])
    assert get_chain_sequence(structure, "L") == ""


def test_get_chain_sequence_none_chain_id_returns_empty():
    structure = _build_structure("H", ["ALA"])
    assert get_chain_sequence(structure, None) == ""


# ---------------------------------------------------------------------------
# _read_fasta_sequences
# ---------------------------------------------------------------------------

def test_read_fasta_sequences_parses_multi_record(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">model.A\nAAAG\nGG\n>model.B\nCCCC\n")
    seqs = _read_fasta_sequences(fasta)
    assert seqs == {"model.A": "AAAGGG", "model.B": "CCCC"}


def test_read_fasta_sequences_empty_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")
    assert _read_fasta_sequences(fasta) == {}


# ---------------------------------------------------------------------------
# correct_json
# ---------------------------------------------------------------------------

def test_correct_json_strips_pdb_suffix(tmp_path):
    annotations_dir = tmp_path / "annotations"
    annotations_dir.mkdir()
    region_json = annotations_dir / "annotations_cdrs.json"
    region_json.write_text(json.dumps({"model_001.pdb": {"cdr": [1, 2]}, "model_002.pdb": {}}))

    result_path = correct_json(tmp_path)

    assert result_path == region_json
    data = json.loads(region_json.read_text())
    assert set(data.keys()) == {"model_001", "model_002"}


def test_correct_json_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        correct_json(tmp_path)


# ---------------------------------------------------------------------------
# hdf5_to_csv
# ---------------------------------------------------------------------------

def test_hdf5_to_csv_converts_predictions(tmp_path):
    hdf5_path = tmp_path / "predictions.hdf5"
    with h5py.File(hdf5_path, "w") as f:
        grp = f.create_group("epoch_0000/pred")
        grp.create_dataset("mol", data=[b"model_001", b"model_002"])
        grp.create_dataset("outputs", data=[0.5, 0.8])

    out_csv = hdf5_to_csv(str(hdf5_path))

    assert out_csv == str(hdf5_path.with_suffix(".csv"))
    df = pd.read_csv(out_csv)
    assert list(df["mol"]) == ["model_001", "model_002"]
    assert list(df["dockq"]) == pytest.approx([0.5, 0.8])


# ---------------------------------------------------------------------------
# ca_coords
# ---------------------------------------------------------------------------

_PDB_LINE = "ATOM  {serial:>5} {name:<4} ALA {chain}{resid:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C\n"


def _write_pdb(path, atoms):
    with open(path, "w") as f:
        for serial, (name, chain, resid, x, y, z) in enumerate(atoms, start=1):
            f.write(
                _PDB_LINE.format(
                    serial=serial, name=name, chain=chain, resid=resid, x=x, y=y, z=z
                )
            )


def test_ca_coords_filters_by_chain_and_atom_name(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    _write_pdb(
        pdb_file,
        [
            ("CA", "H", 1, 1.0, 2.0, 3.0),
            ("CB", "H", 1, 9.0, 9.0, 9.0),  # wrong atom name, skipped
            ("CA", "L", 1, 5.0, 5.0, 5.0),  # wrong chain, skipped
            ("CA", "H", 2, 4.0, 5.0, 6.0),
        ],
    )
    coords = ca_coords(pdb_file, "H")
    assert coords.shape == (2, 3)
    np.testing.assert_allclose(coords, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_ca_coords_respects_max_resid(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    _write_pdb(pdb_file, [("CA", "H", i, float(i), 0.0, 0.0) for i in range(5)])
    coords = ca_coords(pdb_file, "H", max_resid=1)
    assert coords.shape[0] == 2  # residue_index 0 and 1 pass, break when > 1


def test_ca_coords_no_matches_returns_empty_array(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    _write_pdb(pdb_file, [("CA", "H", 1, 1.0, 2.0, 3.0)])
    coords = ca_coords(pdb_file, "Z")
    assert coords.shape == (0, 3)


# ---------------------------------------------------------------------------
# flag_low_contacts
# ---------------------------------------------------------------------------

def test_flag_low_contacts_false_when_missing_chain(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    _write_pdb(pdb_file, [("CA", "H", 1, 0.0, 0.0, 0.0)])
    assert flag_low_contacts(pdb_file, "H", "L") is False


def test_flag_low_contacts_true_when_far_apart(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    # H and L residues placed far apart (>8A) -> zero contacts -> flagged
    atoms = [("CA", "H", i, float(i), 0.0, 0.0) for i in range(1, 4)]
    atoms += [("CA", "L", i, float(i), 100.0, 0.0) for i in range(1, 4)]
    _write_pdb(pdb_file, atoms)
    assert flag_low_contacts(pdb_file, "H", "L", threshold=1) is True


def test_flag_low_contacts_false_when_many_contacts(tmp_path):
    pdb_file = tmp_path / "t.pdb"
    # H and L residues co-located (<8A) -> plenty of contacts -> not flagged
    atoms = [("CA", "H", i, 0.0, 0.0, 0.0) for i in range(1, 20)]
    atoms += [("CA", "L", i, 0.5, 0.0, 0.0) for i in range(1, 20)]
    _write_pdb(pdb_file, atoms)
    assert flag_low_contacts(pdb_file, "H", "L", threshold=15) is False
