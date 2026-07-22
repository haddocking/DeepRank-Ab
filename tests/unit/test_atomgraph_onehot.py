"""Unit test for src/AtomGraph.py's onehot helper.

`onehot` is an instance method but doesn't touch `self`, and `AtomGraph.__init__`
requires a real PDB path (it calls `os.path.basename(pdb)` unconditionally) --
so it's called unbound here rather than via a constructed instance.
"""

import sys
from pathlib import Path

import numpy as np

# AtomGraph.py imports its sibling modules bare (`from tools.edge_orientation
# import ...`), assuming src/ itself is on sys.path -- see CLAUDE.md's
# "Import path inconsistency" note. In the real pipeline this is satisfied as
# a side effect of importing src/GraphGenMP.py first; here we do it directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from src.AtomGraph import AtomGraph


def test_onehot_sets_correct_index():
    result = AtomGraph.onehot(None, idx=2, size=5)
    np.testing.assert_array_equal(result, [0.0, 0.0, 1.0, 0.0, 0.0])


def test_onehot_index_zero():
    result = AtomGraph.onehot(None, idx=0, size=3)
    np.testing.assert_array_equal(result, [1.0, 0.0, 0.0])
