"""Round-trip test for src/Graph.py's nx2h5 / h52nx (NetworkX <-> HDF5).

`type` edge attributes must be bytes (not str) -- nx2h5 stores edge feature
arrays via np.array(...), and h5py has no conversion path for numpy's
unicode ('<U...') dtype, only fixed-width bytes ('S...'). This matches how
the rest of the pipeline (AtomGraph) already writes edge "type" as bytes.
"""

import h5py
import networkx as nx
import numpy as np

from src.Graph import Graph


def _build_sample_graph():
    g = Graph()
    g.name = "test_mol"
    G = nx.DiGraph()
    G.add_node(("A", "1"), x=1.0)
    G.add_node(("A", "2"), x=2.0)
    G.add_node(("B", "1"), x=3.0)
    G.add_edge(("A", "1"), ("B", "1"), type=b"interface", dist=5.0)
    G.add_edge(("A", "1"), ("A", "2"), type=b"internal", dist=1.0)
    g.nx = G
    g.score["dockQ"] = 0.5
    return g


def test_nx2h5_h52nx_roundtrip(tmp_path):
    g = _build_sample_graph()
    hdf5_path = tmp_path / "graph.h5"

    with h5py.File(hdf5_path, "w") as f5:
        g.nx2h5(f5)

    g2 = Graph()
    g2.h52nx(str(hdf5_path), "test_mol")

    assert set(g2.nx.nodes) == set(g.nx.nodes)
    for node in g.nx.nodes:
        assert g2.nx.nodes[node]["x"] == g.nx.nodes[node]["x"]

    assert set(g2.nx.edges) == set(g.nx.edges)
    for u, v in g.nx.edges:
        assert g2.nx.edges[u, v]["type"] == g.nx.edges[u, v]["type"]
        assert g2.nx.edges[u, v]["dist"] == g.nx.edges[u, v]["dist"]

    assert g2.score["dockQ"] == 0.5
    assert g2.score["irmsd"] is None


def test_h52nx_via_existing_group(tmp_path):
    """h52nx also accepts an already-open h5py group (molgrp), not just a path."""
    g = _build_sample_graph()
    hdf5_path = tmp_path / "graph.h5"
    with h5py.File(hdf5_path, "w") as f5:
        g.nx2h5(f5)

    with h5py.File(hdf5_path, "r") as f5:
        g2 = Graph()
        g2.h52nx(None, "test_mol", molgrp=f5["test_mol"])
        assert set(g2.nx.nodes) == set(g.nx.nodes)
        # h5py Group.name is the full in-file path ("/test_mol"), not the bare key
        assert g2.name == "/test_mol"
