import os
import numpy as np
import networkx as nx
import h5py
from pdb2sql import StructureSimilarity
import community
import markov_clustering as mc


class Graph(object):
    """
    Graph-level utilities for protein interface graphs.

    - Compute structural scores (lrmsd, irmsd, fnat, dockQ, bin_class, capri_class)
    - Convert between NetworkX DiGraph and HDF5 layout (nodes, edge types, features)
    - Optional clustering metadata
    """
    def __init__(self):
        self.name = None
        self.nx = None
        self.score = {
            "irmsd": None,
            "lrmsd": None,
            "capri_class": None,
            "fnat": None,
            "dockQ": None,
            "bin_class": None,
        }
        self.clusters = None

    # TODO: function too complex, refactor
    def get_score(self, ref):
        """
        Compute scores against a reference PDB using pdb2sql.StructureSimilarity.
        """
        ref_name = os.path.splitext(os.path.basename(ref))[0]
        sim = StructureSimilarity(self.pdb, ref)

        if os.path.exists(ref_name + ".lzone"):
            self.score["lrmsd"] = sim.compute_lrmsd_fast(method="svd", lzone=ref_name + ".lzone")
            self.score["irmsd"] = sim.compute_irmsd_fast(method="svd", izone=ref_name + ".izone")
        else:
            self.score["lrmsd"] = sim.compute_lrmsd_fast(method="svd")
            self.score["irmsd"] = sim.compute_irmsd_fast(method="svd")

        self.score["fnat"] = sim.compute_fnat_fast()
        self.score["dockQ"] = sim.compute_DockQScore(
            self.score["fnat"], self.score["lrmsd"], self.score["irmsd"]
        )
        self.score["bin_class"] = self.score["irmsd"] < 4.0

        # CAPRI class (lower irmsd -> better class)
        self.score["capri_class"] = 5
        for thr, val in zip([6.0, 4.0, 2.0, 1.0], [4, 3, 2, 1]):
            if self.score["irmsd"] < thr:
                self.score["capri_class"] = val

    def nx2h5(self, f5):
        """
        Write the current NetworkX DiGraph to an HDF5 group layout.
        Preserves directed interface vs. internal edges and their features.
        """
        grp = f5.create_group(self.name)

        # nodes
        data = np.array(list(self.nx.nodes), dtype="S")
        grp.create_dataset("nodes", data=data)

        # node features
        node_feat_grp = grp.create_group("node_data")
        feature_names = list(self.nx.nodes.data())[0][1].keys()
        for feat in feature_names:
            vals = [v for _, v in nx.get_node_attributes(self.nx, feat).items()]
            node_feat_grp.create_dataset(feat, data=np.array(vals))

        # split directed edges by type
        directed_if = [
            (u, v)
            for u, v in self.nx.edges
            if (self.nx.edges[u, v]["type"].decode("utf-8")
                if isinstance(self.nx.edges[u, v]["type"], (bytes, bytearray))
                else self.nx.edges[u, v]["type"]) == "interface"
        ]
        directed_int = [
            (u, v)
            for u, v in self.nx.edges
            if (self.nx.edges[u, v]["type"].decode("utf-8")
                if isinstance(self.nx.edges[u, v]["type"], (bytes, bytearray))
                else self.nx.edges[u, v]["type"]) == "internal"
        ]

        # save textual tuples
        grp.create_dataset("edges", data=np.array(directed_if, dtype="S"))
        grp.create_dataset("internal_edges", data=np.array(directed_int, dtype="S"))

        # integer indices
        node_key = list(self.nx.nodes)
        idx_if = [[node_key.index(u), node_key.index(v)] for u, v in directed_if]
        idx_int = [[node_key.index(u), node_key.index(v)] for u, v in directed_int]
        grp.create_dataset("edge_index", data=np.array(idx_if, dtype=np.int32))
        grp.create_dataset("internal_edge_index", data=np.array(idx_int, dtype=np.int32))

        # edge features (handles 'rij' as stacked vectors)
        edge_feat = list(self.nx.edges(data=True))[0][2].keys()
        edge_grp = grp.create_group("edge_data")
        int_edge_grp = grp.create_group("internal_edge_data")
        for feat in edge_feat:
            vals_if = [self.nx.edges[u, v][feat] for u, v in directed_if]
            vals_int = [self.nx.edges[u, v][feat] for u, v in directed_int]
            if feat == "rij":
                edge_grp.create_dataset("rij", data=np.stack(vals_if, axis=0))
                if vals_int:
                    int_edge_grp.create_dataset("rij", data=np.stack(vals_int, axis=0))
            else:
                edge_grp.create_dataset(feat, data=np.array(vals_if))
                int_edge_grp.create_dataset(feat, data=np.array(vals_int))

        # score
        score_grp = grp.create_group("score")
        for k, v in self.score.items():
            if v is not None:
                score_grp.create_dataset(k, data=v)

        # clustering
        if self.clusters is not None:
            clust_grp = grp.create_group("clustering")
            for method, labels in self.clusters.items():
                method_grp = clust_grp.create_group(method)
                arr = np.asarray(labels)
                if arr.ndim == 1:
                    method_grp.create_dataset("depth_0", data=arr)
                else:
                    for idx in range(arr.shape[1]):
                        method_grp.create_dataset(f"depth_{idx}", data=arr[:, idx])

    def h52nx(self, f5name, mol, molgrp=None):
        """
        Load a NetworkX DiGraph from an HDF5 file/group.
        """
        opened_here = False
        if molgrp is None:
            f5 = h5py.File(f5name, "r")
            molgrp = f5[mol]
            self.name = mol
            self.pdb = mol + ".pdb"
            opened_here = True
        else:
            self.name = molgrp.name
            self.pdb = self.name + ".pdb"

        # directed graph
        self.nx = nx.DiGraph()

        # nodes
        nodes = molgrp["nodes"][()].astype("U").tolist()
        nodes = [tuple(n) for n in nodes]
        node_feats = {k: molgrp[f"node_data/{k}"][()] for k in molgrp["node_data"]}
        for i, n in enumerate(nodes):
            self.nx.add_node(n)
            for k, arr in node_feats.items():
                self.nx.nodes[n][k] = arr[i]

        # interface edges
        edges = [tuple(map(tuple, e)) for e in molgrp["edges"][()].astype("U")]
        edge_feats = {k: molgrp[f"edge_data/{k}"][()] for k in molgrp["edge_data"]}
        for i, (u, v) in enumerate(edges):
            self.nx.add_edge(u, v)
            for k, arr in edge_feats.items():
                self.nx.edges[u, v][k] = arr[i]

        # internal edges
        iedges = [tuple(map(tuple, e)) for e in molgrp["internal_edges"][()].astype("U")]
        if "internal_edge_data" in molgrp:
            i_feats = {k: molgrp[f"internal_edge_data/{k}"][()] for k in molgrp["internal_edge_data"]}
            for i, (u, v) in enumerate(iedges):
                self.nx.add_edge(u, v)
                for k, arr in i_feats.items():
                    self.nx.edges[u, v][k] = arr[i]

        # score
        if "score" in molgrp:
            for k in molgrp["score"]:
                self.score[k] = molgrp[f"score/{k}"][()]

        # clustering
        if "clustering" in molgrp:
            self.clusters = {}
            for method in molgrp["clustering"]:
                gp = molgrp[f"clustering/{method}"]
                depths = [gp[d][()] for d in gp]
                if len(depths) == 1:
                    self.clusters[method] = depths[0]
                else:
                    self.clusters[method] = np.stack(depths, axis=1)

        if opened_here:
            f5.close()

