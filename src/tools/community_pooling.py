#!/usr/bin/env python3
"""
clustering_pooling.py

Utilities to (1) detect communities on residue graphs and (2) pool nodes/edges
by community into a coarsened PyTorch-Geometric graph. Supports Louvain or
Markov Clustering (MCL) for community detection and max/mean aggregation for
features/positions during pooling.
"""

import community
import markov_clustering as mc
import networkx as nx
from scipy.sparse import csr_matrix
import torch

from torch_geometric.utils import scatter
from torch_geometric.nn.pool.pool import pool_edge, pool_batch
from torch_geometric.nn.pool.consecutive import consecutive_cluster
from torch_geometric.data import Batch, Data

import numpy as np
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# -----------------------------------------------------------------------------
# Quick visual sanity check: color nodes by cluster label
# -----------------------------------------------------------------------------
def plot_graph(graph, cluster):
    """Draw a graph with nodes colored by cluster label."""
    import matplotlib.pyplot as plt
    pos = nx.spring_layout(graph, iterations=200)
    nx.draw(graph, pos, node_color=cluster)
    plt.show()


# -----------------------------------------------------------------------------
# Ensure unique cluster labels across graphs in a Batch
# -----------------------------------------------------------------------------
def get_preloaded_cluster(cluster, batch):
    """
    Offset cluster IDs per graph in the batch so labels are globally unique.
    Assumes nodes are ordered by batch.
    """
    nbatch = torch.max(batch) + 1
    for ib in range(1, nbatch):
        cluster[batch == ib] += torch.max(cluster[batch == ib - 1]) + 1
    return cluster


# -----------------------------------------------------------------------------
# Community detection
# -----------------------------------------------------------------------------
def community_detection_per_batch(
    edge_index, batch, num_nodes, edge_attr=None, method="mcl"
):
    """Detect clusters independently per graph in a Batch.

    Builds a NetworkX graph for each batch slice, applies the chosen method,
    and returns a single tensor of cluster labels with globally unique IDs.
    """
    # Build global NetworkX graph once
    g = nx.Graph()
    g.add_nodes_from(range(num_nodes))
    for idx, (i, j) in enumerate(edge_index.transpose(0, 1).tolist()):
        if edge_attr is None:
            g.add_edge(i, j)
        else:
            g.add_edge(i, j, weight=edge_attr[idx])

    num_batch = int(batch.max()) + 1
    all_idx = list(range(num_nodes))
    clusters_out = []
    ncluster = 0

    for b in range(num_batch):
        # nodes for batch b, then deterministic order
        nodes_b = [n for n in all_idx if batch[n] == b]
        subg = g.subgraph(nodes_b)
        nodes_sorted = sorted(subg.nodes())

        if method == "louvain":
            part = community.best_partition(subg)
            labels = [part[n] for n in nodes_sorted]
            for lbl in labels:
                clusters_out.append(lbl + ncluster)
            ncluster = max(clusters_out) + 1

        elif method == "mcl":
            matrix = nx.to_scipy_sparse_matrix(subg, nodelist=nodes_sorted)
            result = mc.run_mcl(matrix)
            mcl_clusters = mc.get_clusters(result)
            arr = np.zeros(len(nodes_sorted), dtype=int)
            for cid, comp in enumerate(mcl_clusters):
                for idx in comp:
                    arr[idx] = cid + ncluster
            clusters_out.extend(arr.tolist())
            ncluster = clusters_out[-1] + 1

        else:
            raise ValueError(f"clustering method {method} not supported")

    device = edge_index.device
    return torch.tensor(clusters_out, device=device)


def community_detection(edge_index, num_nodes, edge_attr=None, method="mcl"):
    """Detect clusters on a single graph with deterministic node ordering."""
    g = nx.Graph()
    g.add_nodes_from(range(num_nodes))
    for idx, (i, j) in enumerate(edge_index.transpose(0, 1).tolist()):
        if edge_attr is None:
            g.add_edge(i, j)
        else:
            g.add_edge(i, j, weight=edge_attr[idx])

    device = edge_index.device
    nodes_sorted = sorted(g.nodes())

    if method == "louvain":
        part = community.best_partition(g)
        labels = [part[n] for n in nodes_sorted]
        return torch.tensor(labels, device=device)

    elif method == "mcl":
        matrix = csr_matrix(nx.adjacency_matrix(g, nodelist=nodes_sorted))
        result = mc.run_mcl(matrix)
        mcl_clusters = mc.get_clusters(result)
        arr = np.zeros(num_nodes, dtype=int)
        for cid, comp in enumerate(mcl_clusters):
            for idx in comp:
                arr[idx] = cid
        return torch.tensor(arr, device=device)

    else:
        raise ValueError(f"clustering method {method} not supported")


# -----------------------------------------------------------------------------
# Pooling by community
# -----------------------------------------------------------------------------
def community_pooling(cluster, data):
    """Pool nodes/edges by community into a coarsened PyG graph.

    Args:
        cluster (LongTensor[N]): cluster id per node.
        data (Data or Batch): graph with
            - x: (s, v) tuple or a single tensor
            - edge_index, edge_attr: (es, ev)
            - optional: pos, pos2D, batch, cluster0/cluster1
            - optional: internal_edge_index, internal_edge_attr (tuple or tensor)
    Returns:
        Data or Batch: coarsened graph with pooled features and edges.
    """
    # 1) normalize cluster ids and align to device
    cluster, perm = consecutive_cluster(cluster)
    x0 = data.x[0] if isinstance(data.x, tuple) else data.x
    cluster = cluster.to(x0.device)

    # 2) pool node features: max over scalars, mean over vectors
    if isinstance(data.x, tuple):
        s, v = data.x
        new_s = scatter(s, cluster, dim=0, reduce="max")
        new_v = scatter(v, cluster, dim=0, reduce="mean") if v is not None else None
        new_x = (new_s, new_v)
    else:
        new_s = scatter(data.x, cluster, dim=0, reduce="max")
        new_x = new_s

    # 3) pool external edges (separate scalar/vector parts)
    es, ev = data.edge_attr
    ei = data.edge_index
    new_ei, new_es = pool_edge(cluster, ei, es)
    _,       new_ev = pool_edge(cluster, ei, ev)
    new_edge_attr = (new_es, new_ev)

    # 4) pool internal edges if present
    if hasattr(data, "internal_edge_index") and data.internal_edge_index is not None:
        int_ei = data.internal_edge_index
        orig_int_attr = data.internal_edge_attr
        if isinstance(orig_int_attr, tuple):
            int_es, int_ev = orig_int_attr
        else:
            int_es, int_ev = orig_int_attr, None

        new_int_ei, new_int_es = pool_edge(cluster, int_ei, int_es)
        if int_ev is not None:
            _, new_int_ev = pool_edge(cluster, int_ei, int_ev)
        else:
            new_int_ev = None
        new_internal_edge_index = new_int_ei
        new_internal_edge_attr  = (new_int_es, new_int_ev)
    else:
        new_internal_edge_index = None
        new_internal_edge_attr  = None

    # 5) pool positions if present
    new_pos   = scatter(data.pos,   cluster, dim=0, reduce="mean") if hasattr(data, "pos")   else None
    new_pos2D = scatter(data.pos2D, cluster, dim=0, reduce="mean") if hasattr(data, "pos2D") else None

    # 6) pool batch if present
    if hasattr(data, "batch") and data.batch is not None:
        new_batch = pool_batch(perm, data.batch.to(cluster.device))
    else:
        new_batch = None

    # 7) build output object
    if isinstance(data, Batch):
        out = Batch(
            batch=new_batch,
            x=new_x,
            edge_index=new_ei,
            edge_attr=new_edge_attr,
        )
    else:
        out = Data(
            x=new_x,
            edge_index=new_ei,
            edge_attr=new_edge_attr,
        )

    # 8) attach optional attributes
    if new_internal_edge_index is not None:
        out.internal_edge_index = new_internal_edge_index
        out.internal_edge_attr  = new_internal_edge_attr
    if new_pos is not None:
        out.pos = new_pos
    if new_pos2D is not None:
        out.pos2D = new_pos2D
    if new_batch is not None:
        out.batch = new_batch
    if hasattr(data, "cluster0"):
        out.cluster0 = data.cluster0
        out.cluster1 = data.cluster1

    return out
