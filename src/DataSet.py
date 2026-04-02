import sys
import os
import torch
import numpy as np
from torch_geometric.data.dataset import Dataset
from torch_geometric.data.data import Data
from tqdm import tqdm
import h5py
import copy
from sklearn.preprocessing import RobustScaler

from src.tools.community_pooling import community_detection, community_pooling

"""
Load graphs from one or more HDF5 files and produce PyG `Data` objects.

- Reads precomputed node/edge features per complex.
- Supports on-the-fly indexing, filtering, and simple train/val splits.
- Optionally precomputes and caches hierarchical clustering (depth 0/1).
"""

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# ---------------------------
# Utilities
# ---------------------------

def DivideDataSet(dataset, percent=[0.8, 0.2], shuffle=True):
    """
    Split the dataset into two parts by cloning indices.
    """
    size = dataset.__len__()
    index = np.arange(size)
    if shuffle:
        np.random.shuffle(index)

    size1 = int(percent[0] * size)
    index1, index2 = index[:size1], index[size1:]

    ds1 = copy.deepcopy(dataset)
    ds1.index_complexes = [dataset.index_complexes[i] for i in index1]

    ds2 = copy.deepcopy(dataset)
    ds2.index_complexes = [dataset.index_complexes[i] for i in index2]

    return ds1, ds2


def PreCluster(dataset, method):
    """
    Precompute and cache node clustering (mcl or louvain) at two depths.
    """
    for fname, mol in tqdm(dataset.index_complexes):
        data = dataset.load_one_graph(fname, mol)

        if data is None:
            f5 = h5py.File(fname, "a")
            try:
                print(f"deleting {mol}")
                del f5[mol]
            except Exception:
                print(f"{mol} not found")
            f5.close()
            continue

        f5 = h5py.File(fname, "a")
        grp = f5[mol]
        clust_grp = grp.require_group("clustering")

        if method.lower() not in clust_grp:
            method_grp = clust_grp.create_group(method.lower())

            cluster = community_detection(data.internal_edge_index, data.num_nodes, method=method)
            method_grp.create_dataset("depth_0", data=cluster)

            data = community_pooling(cluster, data)
            cluster = community_detection(data.internal_edge_index, data.num_nodes, method=method)
            method_grp.create_dataset("depth_1", data=cluster)

            f5.close()
        else:
            f5.close()
            continue


# ---------------------------
# Dataset
# ---------------------------

class HDF5DataSet(Dataset):
    def __init__(
        self,
        root="./",
        name=None,
        database=None,
        transform=None,
        pre_transform=None,
        dict_filter=None,
        target=None,
        tqdm=True,
        index=None,
        node_feature="all",
        edge_feature="all",
        clustering_method="mcl",
        dist_transform=lambda x: np.tanh(-x / 2 + 2) + 1,
        rij_transform=lambda arr: arr / (np.linalg.norm(arr, axis=1, keepdims=True) + 1e-6),
    ):
        """
        Initialize dataset backed by one or more HDF5 files.
        """
        super().__init__(root, transform, pre_transform)

        # accept a single file or a list
        self.database = database if isinstance(database, list) else [database]
        self.name = name
        self.target = target
        self.dict_filter = dict_filter
        self.tqdm = tqdm
        self.index = index

        self.node_feature = node_feature
        self.edge_feature = edge_feature


        self.dist_transform = dist_transform
        self.rij_transform = rij_transform
        self.clustering_method = clustering_method

        # file integrity and feature schema
        self.check_hdf5_files()
        self.check_node_feature()
        self.check_edge_feature()

        # build master index of (file, group)
        self.create_index_molecules()

        # fit a RobustScaler on edge feature "elec" if present
        elec_vals = []
        for fname, mol in self.index_complexes:
            with h5py.File(fname, "r") as f5:
                grp = f5[mol]
                if "elec" in grp["edge_data"]:
                    arr = grp["edge_data"]["elec"][()]
                    elec_vals.append(arr.reshape(-1, 1))
        if elec_vals:
            all_elec = np.vstack(elec_vals)
            self.elec_scaler = RobustScaler(
                quantile_range=(25, 75),
                with_centering=True,
                with_scaling=True,
            ).fit(all_elec)
        else:
            self.elec_scaler = None

    def len(self):
        """
        Number of complexes.
        """
        return len(self.index_complexes)

    def _download(self):
        pass

    def _process(self):
        pass

    def get(self, index):
        """
        Load one graph by index.
        """
        fname, mol = self.index_complexes[index]
        return self.load_one_graph(fname, mol)

    # ---------------------------
    # Schema checks
    # ---------------------------

    def check_hdf5_files(self):
        """
        Basic integrity check: drop empty or unreadable files.
        """
        remove_file = []
        for fname in self.database:
            try:
                f = h5py.File(fname, "r")
                if len(list(f.keys())) == 0:
                    print(f"empty file {fname}")
                    remove_file.append(fname)
                f.close()
            except Exception as e:
                print(e)
                print(f"corrupted file {fname}")
                remove_file.append(fname)
        for name in remove_file:
            self.database.remove(name)

    def check_node_feature(self):
        """
        Collect available node features from the first file; validate selection.
        """
        f = h5py.File(self.database[0], "r")
        mol_key = list(f.keys())[0]
        self.available_node_feature = list(f[mol_key + "/node_data/"].keys())
        f.close()

        if self.node_feature == "all":
            self.node_feature = self.available_node_feature
        else:
            for feat in self.node_feature:
                if feat not in self.available_node_feature:
                    print(feat, "node feature not found in", self.database[0])
                    print("available node features:")
                    print("\n".join(self.available_node_feature))
                    exit()

    def check_edge_feature(self):
        """
        Do not pre-type features.

        Scalar:
          - any key in /edge_data can be requested (use 'all' for all keys).

        Vector:
          - 'all' means keys whose last dimension looks like 3 (heuristic).
          - explicit requests are only validated for existence here.
        """
        with h5py.File(self.database[0], "r") as f:
            mol_key = list(f.keys())[0]
            edge_grp = f[f"{mol_key}/edge_data/"]
            all_keys = [k for k in edge_grp.keys()]

            if self.edge_feature == "all":
                self.edge_feature = list(all_keys)
            else:
                for feat in self.edge_feature or []:
                    if feat not in all_keys:
                        raise ValueError(f"scalar edge feature {feat!r} not found; available: {all_keys}")

            vec_keys_guess = []
            for k in all_keys:
                try:
                    shp = edge_grp[k].shape
                except Exception:
                    continue
                if isinstance(shp, tuple) and len(shp) >= 2 and shp[-1] == 3:
                    vec_keys_guess.append(k)



        self.available_edge_feature = self.edge_feature

    # ---------------------------
    # Loader
    # ---------------------------

    def load_one_graph(self, fname, mol):
        import h5py, numpy as np, torch
        from torch_geometric.data import Data

        # convert any /edge_data array to (M, K) for scalar use
        def _as_scalar_matrix(arr: np.ndarray, edge_index: torch.Tensor, pos: torch.Tensor) -> np.ndarray:
            """
            Normalize shapes for scalar edge features.
            (M,) -> (M,1); (M,K) -> (M,K); (M,3) -> (M,3); (M,F,3) -> (M, 3F).
            For vector-like inputs, disallow (E,2,3); (E,3) is allowed with E or 2E.
            """
            M = edge_index.size(1)

            if isinstance(arr, np.ndarray) and arr.ndim >= 2 and arr.shape[-1] == 3:
                if arr.ndim == 3 and arr.shape[1] == 2 and arr.shape[2] == 3:
                    raise RuntimeError("found (E,2,3) but pair alignment is disabled")

                if arr.ndim == 3 and arr.shape[0] == M:
                    return arr.reshape(M, -1)

                if arr.ndim == 2 and arr.shape[1] == 3:
                    E = arr.shape[0]
                    if M == 2 * E:
                        return np.repeat(arr, 2, axis=0)
                    elif M == E:
                        return arr
                    else:
                        raise RuntimeError(f"scalar-usage: (E,3) cannot match M={M} (E={E})")

                raise RuntimeError(f"scalar-usage: unexpected vector-like shape {arr.shape}")

            if arr.ndim == 1:
                return arr.reshape(-1, 1)
            if arr.ndim == 2:
                return arr

            raise RuntimeError(f"scalar-usage: unexpected shape {arr.shape}")

        with h5py.File(fname, "r") as f5:
            if mol not in f5:
                return None
            grp = f5[mol]

            # nodes
            node_scalars = []
            for feat in self.node_feature:
                if feat == 'embedding':
                    #repalce the embedding with mean along the embedding dimension
                    arr = grp[f"node_data/{feat}"][()].mean(axis=1)
                else:
                    arr = grp[f"node_data/{feat}"][()]
                node_scalars.append(arr.reshape(-1, 1) if arr.ndim == 1 else arr)
            node_s = torch.tensor(np.hstack(node_scalars), dtype=torch.float32)

            pos = torch.tensor(grp["node_data/pos"][()], dtype=torch.float32)
            node_v = pos.unsqueeze(1)

            # external (directed)
            edge_index = torch.tensor(grp["edge_index"][()].T, dtype=torch.long).contiguous()
            M = edge_index.size(1)

            # scalar external edges
            edge_scalars = []
            edge_data_keys = set(grp["edge_data"].keys()) if "edge_data" in grp else set()
            for feat in self.edge_feature:
                if feat not in edge_data_keys:
                    continue
                arr = grp[f"edge_data/{feat}"][()]
                if not np.issubdtype(arr.dtype, np.number):
                    continue

                try:
                    arr = _as_scalar_matrix(arr, edge_index, pos)
                except RuntimeError as e:
                    raise RuntimeError(f"{e} in {os.path.basename(fname)}::{mol} for scalar feature {feat!r}")

                if feat == "dist":
                    arr = self.dist_transform(arr)
                elif feat == "elec" and self.elec_scaler is not None:
                    arr = self.elec_scaler.transform(arr)
                elif feat == "rij":
                    arr = self._normalize_rij_scalar(arr)

                edge_scalars.append(arr)

            edge_s = torch.tensor(np.hstack(edge_scalars), dtype=torch.float32) if edge_scalars else None
            S_ext = edge_s.shape[1] if edge_s is not None else 0

            # vector external edges
            vecs = []
            edge_data_keys = set(grp["edge_data"].keys()) if "edge_data" in grp else set()

            edge_v = torch.tensor(np.concatenate(vecs, axis=1), dtype=torch.float32) if vecs else torch.zeros((M, 0, 3), dtype=torch.float32)
            F_ext = edge_v.shape[1]

            # internal edges
            arr_ie = grp["internal_edge_index"][()] if "internal_edge_index" in grp else None
            if arr_ie is None or (isinstance(arr_ie, np.ndarray) and arr_ie.ndim == 1):
                internal_edge_index = torch.empty((2, 0), dtype=torch.long)
            else:
                internal_edge_index = torch.tensor(arr_ie.T, dtype=torch.long).contiguous()
            Mi = internal_edge_index.size(1)

            # scalar internal edges
            int_scalars = []
            int_edge_scalar_keys = set(grp["internal_edge_data"].keys()) if "internal_edge_data" in grp else set()
            for feat in self.edge_feature:
                if feat not in int_edge_scalar_keys:
                    continue
                arr_i = grp[f"internal_edge_data/{feat}"][()]
                if not np.issubdtype(arr_i.dtype, np.number):
                    continue

                try:
                    arr_i = _as_scalar_matrix(arr_i, internal_edge_index, pos)
                except RuntimeError as e:
                    raise RuntimeError(f"{e} in {os.path.basename(fname)}::{mol} (internal) for scalar feature {feat!r}")

                if feat == "dist":
                    arr_i = self.dist_transform(arr_i)
                elif feat == "elec" and self.elec_scaler is not None and arr_i.shape[0] > 0:
                    arr_i = self.elec_scaler.transform(arr_i)
                elif feat == "rij":
                    arr_i = self._normalize_rij_scalar(arr_i)

                int_scalars.append(arr_i)

            internal_edge_attr = torch.tensor(np.hstack(int_scalars), dtype=torch.float32) if int_scalars else torch.zeros((Mi, 0), dtype=torch.float32)

            # match internal scalar width to external scalar width
            S_int = internal_edge_attr.shape[1]
            if S_ext != S_int:
                if S_int < S_ext:
                    pad = torch.zeros((Mi, S_ext - S_int), dtype=torch.float32)
                    internal_edge_attr = torch.cat([internal_edge_attr, pad], dim=1)
                else:
                    internal_edge_attr = internal_edge_attr[:, :S_ext]

            # vector internal edges
            int_vecs = []
            int_edge_vec_keys = set(grp["internal_edge_data"].keys()) if "internal_edge_data" in grp else set()
            for vf in []:
                if vf in int_edge_vec_keys:
                    arr_i = grp[f"internal_edge_data/{vf}"][()]
                    if arr_i.ndim == 1 and arr_i.size == 0:
                        int_vecs.append(np.zeros((Mi, 1, 3), dtype=np.float32))
                        continue

                    if arr_i.ndim == 3 and arr_i.shape[1] == 2 and arr_i.shape[2] == 3:
                        raise RuntimeError(f"internal vec-edge {vf}: found (E,2,3) but pair alignment is disabled in {os.path.basename(fname)}::{mol}")

                    if vf == "rij":
                        arr_i = self.rij_transform(arr_i)

                    if arr_i.ndim == 2 and arr_i.shape[1] == 3:
                        arr_i = arr_i[:, None, :]
                    elif arr_i.ndim == 3 and arr_i.shape[-1] == 3:
                        pass
                    else:
                        raise RuntimeError(f"unexpected internal vec-edge shape {arr_i.shape} for {vf} in {mol}")

                    if arr_i.shape[0] != Mi:
                        raise RuntimeError(f"internal vec-edge {vf}: length {arr_i.shape[0]} != internal_edge_index {Mi} in {os.path.basename(fname)}::{mol}")

                    int_vecs.append(arr_i)
                else:
                    int_vecs.append(np.zeros((Mi, 1, 3), dtype=np.float32))

            internal_edge_vec = torch.tensor(np.concatenate(int_vecs, axis=1), dtype=torch.float32) if int_vecs else torch.zeros((Mi, 0, 3), dtype=torch.float32)

            # match internal vector width to external vector width
            F_int = internal_edge_vec.shape[1]
            if F_ext != F_int:
                if F_int < F_ext:
                    pad = torch.zeros((Mi, F_ext - F_int, 3), dtype=torch.float32)
                    internal_edge_vec = torch.cat([internal_edge_vec, pad], dim=1)
                else:
                    internal_edge_vec = internal_edge_vec[:, :F_ext, :]

            # quick checks
            assert edge_v.shape[0] == edge_index.size(1)
            assert edge_v.shape[-1] == 3
            assert internal_edge_vec.shape[-1] == 3
            assert internal_edge_attr.shape[1] == S_ext

            # target
            y = None
            if self.target and "score" in grp and self.target in grp["score"]:
                y = torch.tensor([grp[f"score/{self.target}"][()]], dtype=torch.float32)

            # assemble
            data = Data(
                x=(node_s, node_v),
                edge_index=edge_index,
                edge_attr=(edge_s, edge_v),
                y=y,
                pos=pos,
            )
            data.internal_edge_index = internal_edge_index
            data.internal_edge_attr = (internal_edge_attr, internal_edge_vec)
            data.mol = mol

            if "clustering" in grp and self.clustering_method in grp["clustering"]:
                clgrp = grp[f"clustering/{self.clustering_method}"]
                if "depth_0" in clgrp and "depth_1" in clgrp:
                    data.cluster0 = torch.tensor(clgrp["depth_0"][()], dtype=torch.long)
                    data.cluster1 = torch.tensor(clgrp["depth_1"][()], dtype=torch.long)

            return data

    # ---------------------------
    # Helpers
    # ---------------------------

    def _normalize_rij_scalar(self, arr: np.ndarray) -> np.ndarray:
        """
        Apply `rij_transform` and return 2D shape for scalar use:
        (E, 3) -> (E, 3); (E, F, 3) -> (E, 3F).
        """
        if arr.ndim == 2 and arr.shape[1] == 3:
            return self.rij_transform(arr)
        if arr.ndim == 3 and arr.shape[-1] == 3:
            E, F, _ = arr.shape
            flat = arr.reshape(E * F, 3)
            flat = self.rij_transform(flat).reshape(E, F, 3)
            return flat.reshape(E, -1)
        return arr

    # ---------------------------
    # Indexing
    # ---------------------------

    def create_index_molecules(self):
        """
        Build the (file, group) index for all complexes, with optional filtering.
        """
        self.index_complexes = []

        desc = "{:25s}".format("train dataset")
        data_iter = tqdm(self.database, desc=desc, file=sys.stdout) if self.tqdm else self.database
        if not self.tqdm:
            sys.stdout.flush()

        for fdata in data_iter:
            if self.tqdm:
                data_iter.set_postfix(mol=os.path.basename(fdata))
            try:
                fh5 = h5py.File(fdata, "r")
                if self.index is None:
                    mol_names = list(fh5.keys())
                else:
                    mol_names = [list(fh5.keys())[i] for i in self.index]

                for k in mol_names:
                    if self.filter(fh5[k]):
                        self.index_complexes.append((fdata, k))
                fh5.close()
            except Exception as inst:
                print(f"ignore file {fdata}")
                print(inst)

        self.ntrain = len(self.index_complexes)
        self.index_train = list(range(self.ntrain))
        self.ntot = len(self.index_complexes)

    def filter(self, molgrp):
        """
        Optional filter based on `self.dict_filter`, e.g., {'dockq': '>0.5'}.
        """
        if self.dict_filter is None:
            return True

        for cond_name, cond_vals in self.dict_filter.items():
            try:
                val = molgrp["score"][cond_name][()]
            except KeyError:
                print(f"filter {cond_name} not found")
                print("available keys:")
                for k in molgrp["score"].keys():
                    print(k)
                return False

            if isinstance(cond_vals, str):
                expr = cond_vals
                for o in [">", "<", "=="]:
                    expr = expr.replace(o, f"val{o}")
                if not eval(expr):
                    return False
            else:
                raise ValueError("conditions not supported", cond_vals)

        return True
