#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EGNN-based model with dual (external/internal) branches.

Key points
- Edge vector features are used directly by flattening (E, F, 3) -> (E, F*3) and
  concatenating to scalar edge features.
- External API remains the same: `edge_in_v` denotes F (number of vector features per edge).
- Comments are high-level; brief inline notes appear only where operations are non-obvious.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_scatter import scatter_add
from torch_geometric.nn import GlobalAttention, max_pool_x
from torch_geometric.utils import dropout_adj
from tools.community_pooling import get_preloaded_cluster, community_pooling

#ignore warnings
import warnings
warnings.filterwarnings(
    "ignore",
    message="'nn.glob.GlobalAttention' is deprecated",
    category=UserWarning
)



# ---- utilities ----

class GaussianRBF(nn.Module):
    """Radial basis functions over distances to stabilize learning on |r_ij|."""
    def __init__(self, n_rbf=16, r_min=0.0, r_max=30.0):
        super().__init__()
        self.centers = nn.Parameter(torch.linspace(r_min, r_max, n_rbf), requires_grad=False)
        self.widths  = nn.Parameter(torch.tensor((r_max - r_min) / n_rbf), requires_grad=False)

    def forward(self, d: torch.Tensor) -> torch.Tensor:
        x = (d.unsqueeze(-1) - self.centers) / (self.widths + 1e-8)
        return torch.exp(-x**2)  # (E, n_rbf)


def _split_attr_maybe_tuple(attr):
    """Return (scalars, vectors); handle tuple/list/single tensor/None."""
    if isinstance(attr, (tuple, list)):
        if len(attr) == 2:
            return attr[0], attr[1]
        elif len(attr) == 1:
            return attr[0], None
        else:
            return None, None
    return attr, None


# ---- EGNN blocks ----

class EGNNBlockBalanced(nn.Module):
    """EGNN block using flattened edge vectors (no vector stats)."""

    def __init__(self, s_dim, e_s_dim, e_v_flat_dim=0, n_rbf=18, dropout=0.08):
        """
        e_v_flat_dim should be edge_in_v * 3 when vector edge features are present,
        else 0.
        """
        super().__init__()
        self.rbf = GaussianRBF(n_rbf=n_rbf)
        self.e_v_flat_dim = int(e_v_flat_dim)

        # message input: h_src, h_dst, rbf(d_ij), edge scalars, flattened edge vectors
        in_e = 2*s_dim + n_rbf + e_s_dim + self.e_v_flat_dim

        self.phi_e = nn.Sequential(
            nn.Linear(in_e, 320), nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(320, 160), nn.SiLU(),
            nn.Dropout(dropout * 0.5),
            nn.Linear(160, 128), nn.SiLU()
        )

        self.gate = nn.Sequential(
            nn.Linear(128, 64), nn.ReLU(),
            nn.Linear(64, 1), nn.Sigmoid()
        )

        self.phi_h = nn.Sequential(
            nn.Linear(128, s_dim), nn.SiLU()
        )

        self.phi_x = nn.Sequential(
            nn.Linear(128, 32), nn.SiLU(),
            nn.Linear(32, 1)
        )

        self.ln = nn.LayerNorm(s_dim)
        self.alpha = nn.Parameter(torch.tensor(0.8))

    def forward(self, h, x, edge_index, e_s=None, e_v=None, x_scale=0.08):
        if edge_index is None or edge_index.numel() == 0:
            return h, x

        src, dst = edge_index
        r_ij = x[src] - x[dst]
        d_ij = torch.norm(r_ij, dim=-1)
        rbf = self.rbf(d_ij)

        if e_s is None:
            e_s = torch.zeros((r_ij.size(0), 0), device=h.device, dtype=h.dtype)

        # flatten vector edge features (E, F, 3) -> (E, F*3); keep empty tensor if absent
        if self.e_v_flat_dim == 0:
            ev_flat = torch.zeros((r_ij.size(0), 0), device=h.device, dtype=h.dtype)
        else:
            if (e_v is None) or (not isinstance(e_v, torch.Tensor)) or (e_v.numel() == 0):
                ev_flat = torch.zeros((r_ij.size(0), self.e_v_flat_dim), device=h.device, dtype=h.dtype)
            else:
                ev_flat = e_v.reshape(e_v.size(0), -1)
                
        m_in = torch.cat([h[src], h[dst], rbf, e_s, ev_flat], dim=1)

        m = self.phi_e(m_in)
        a = self.gate(m)
        m = m * a

        dh = self.phi_h(m)
        dh_agg = scatter_add(dh, dst, dim=0, dim_size=h.size(0))
        h_new = self.ln(h + torch.sigmoid(self.alpha) * dh_agg)

        # coordinate update along r_ij with learned scalar coefficient
        coeff = self.phi_x(m) * x_scale
        dx = r_ij * coeff
        dx_agg = scatter_add(dx, dst, dim=0, dim_size=x.size(0))
        x_new = x + dx_agg

        return h_new, x_new


# ---- Main Network ----

class egnn(nn.Module):
    """
    Uses flattened edge vector features when edge_in_v > 0.
    """
    def __init__(self, node_in_s, node_in_v, edge_in_s, edge_in_v, out_dim=1):
        """
        node_in_v/edge_in_v can be 0 when vector features are not used.
        edge_in_v = F (number of vector features per edge). Internally uses F*3 as flat input.
        """
        super().__init__()

        self.use_edge_vec = bool(edge_in_v and edge_in_v > 0)
        self._did_info_print = False

        # balanced hyperparameters
        self.n_layers = 7
        self.s_hidden = 180
        self.n_rbf = 18

        # dropout settings
        self.dropedge_p = 0.15
        self.dropnode_p = 0.12
        self.layer_dropout = 0.08
        self.head_dropout = 0.25

        self.x_scale = 0.08

        # node scalar encoder (node vectors are not used in this architecture)
        self.s_in = nn.Sequential(
            nn.LayerNorm(node_in_s),
            nn.Linear(node_in_s, self.s_hidden // 2),
            nn.SiLU(),
            nn.Dropout(0.05),
            nn.Linear(self.s_hidden // 2, self.s_hidden),
            nn.SiLU()
        )

        # flattened edge vector feature dimension
        e_v_flat_dim = (edge_in_v * 3) if self.use_edge_vec else 0

        def make_stack():
            return nn.ModuleList([
                EGNNBlockBalanced(
                    self.s_hidden,
                    e_s_dim=edge_in_s,
                    e_v_flat_dim=e_v_flat_dim,
                    n_rbf=self.n_rbf,
                    dropout=self.layer_dropout
                )
                for _ in range(self.n_layers)
            ])

        self.blocks_ext = make_stack()
        self.blocks_int = make_stack()

        pooled_dim = self.s_hidden + 1

        def make_pool():
            return GlobalAttention(
                gate_nn=nn.Sequential(
                    nn.Linear(pooled_dim, 96), nn.ReLU(),
                    nn.Dropout(0.1),
                    nn.Linear(96, 1)
                ),
                nn=nn.Sequential(
                    nn.Linear(pooled_dim, 192), nn.ReLU(),
                    nn.Dropout(0.1),
                    nn.Linear(192, 96)
                )
            )

        self.pool_ext = make_pool()
        self.pool_int = make_pool()

        self.fc1 = nn.Linear(192, 128)   # 96*2 -> 128
        self.fc_res = nn.Linear(192, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, out_dim)

        self.apply(self._init_weights)

    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight, gain=0.4)
            if module.bias is not None:
                nn.init.constant_(module.bias, 0)

    @staticmethod
    def _dropnode(h, x, p):
        if p <= 0 or not h.requires_grad:
            return h, x
        mask = (torch.rand(h.size(0), device=h.device) > p).float().unsqueeze(1)
        return h * mask, x * mask

    @staticmethod
    def _dropedge(ei, e_s, e_v, p):
        """
        Works with or without vector edge features. When e_v is (E, F, 3),
        flatten to (E, F*3) to apply dropout_adj, then restore shape.
        """
        if ei is None or p <= 0 or ei.numel() == 0:
            return ei, e_s, e_v
        es_dim = e_s.size(1) if e_s is not None else 0
        ev_dim = e_v.size(1) if (e_v is not None and isinstance(e_v, torch.Tensor) and e_v.numel() > 0) else 0

        if e_v is None:
            if e_s is None:
                return ei, e_s, e_v
            e_v = torch.zeros((e_s.size(0), 0, 3), device=e_s.device, dtype=e_s.dtype)
        if e_s is None:
            e_s = torch.zeros((e_v.size(0), 0), device=e_v.device, dtype=e_v.dtype)

        ev_flat = e_v.view(e_v.size(0), -1)  # (E, F*3) or (E, 0)
        ei2, comb = dropout_adj(ei, torch.cat([e_s, ev_flat], dim=1), p=p, force_undirected=True, training=True)
        e_s2 = comb[:, :es_dim]
        if ev_dim > 0:
            e_v2 = comb[:, es_dim:].view(-1, ev_dim, 3)
        else:
            e_v2 = torch.zeros((comb.size(0), 0, 3), device=comb.device, dtype=comb.dtype)
        return ei2, e_s2, e_v2

    @staticmethod
    def _hier_pool(data_like, h, x):
        # package (scalars, vectors) for community_pooling; vectors as 1-channel
        v = x.unsqueeze(1)
        data_like.x = (h, v)

        # first level pooling
        cl0 = get_preloaded_cluster(data_like.cluster0, data_like.batch)
        data_like = community_pooling(cl0, data_like)
        h, v = data_like.x
        x = v.squeeze(1)

        # second level pooling; concatenate norm of vector as an extra scalar
        cl1 = get_preloaded_cluster(data_like.cluster1, data_like.batch)
        vnorm = torch.norm(v, dim=-1).squeeze(1)
        z = torch.cat([h, vnorm.unsqueeze(1)], dim=1)
        x_pool, b_pool = max_pool_x(cl1, z, data_like.batch)
        return x_pool, b_pool

    def forward(self, data):
        # ----- inputs -----
        s_in, _v_in = _split_attr_maybe_tuple(getattr(data, "x", None))
        if s_in is None:
            raise ValueError("data.x missing or malformed: expected scalar tensor or (scalars, vectors).")

        e_s_ext, e_v_ext = _split_attr_maybe_tuple(getattr(data, "edge_attr", None))
        ei_ext = getattr(data, "edge_index", None)

        iei = getattr(data, "internal_edge_index", None)
        e_s_int, e_v_int = _split_attr_maybe_tuple(getattr(data, "internal_edge_attr", None))

        if not hasattr(data, "pos") or getattr(data, "pos") is None:
            raise ValueError("data.pos missing: coordinates are required.")
        x0 = data.pos.clone()

        if not self._did_info_print:
            def _present(t):
                return (t is not None) and isinstance(t, torch.Tensor) and (t.numel() > 0)
            ext_state = "present" if _present(e_v_ext) else "absent"
            int_state = "present" if _present(e_v_int) else "absent"
            self._did_info_print = True

        # ----- pipeline -----
        h0 = self.s_in(s_in)
        x = x0

        if self.training:
            h0, x = self._dropnode(h0, x, self.dropnode_p)

        # ===== external branch =====
        ei_e, es_e, ev_e = ei_ext, e_s_ext, (e_v_ext if self.use_edge_vec else None)
        if self.training:
            ei_e, es_e, ev_e = self._dropedge(ei_e, es_e, ev_e, self.dropedge_p)

        h_e, x_e = h0, x
        h_skip = h_e
        for i, blk in enumerate(self.blocks_ext):
            # pass ev_e directly; the block flattens it
            h_e, x_e = blk(h_e, x_e, ei_e, es_e, ev_e, x_scale=self.x_scale)
            if (i + 1) % 2 == 0 and i > 0:
                h_e = h_e + 0.2 * h_skip
                h_skip = h_e

        data_e = data.clone()
        data_e.edge_index = ei_e
        data_e.edge_attr = (es_e, ev_e)
        data_e.internal_edge_index = iei
        data_e.internal_edge_attr = (e_s_int, (e_v_int if self.use_edge_vec else None))
        xpool_e, bpool_e = self._hier_pool(data_e, h_e, x_e)
        g_e = self.pool_ext(xpool_e, bpool_e)

        # ===== internal branch =====
        ei_i, es_i, ev_i = iei, e_s_int, (e_v_int if self.use_edge_vec else None)
        if ei_i is not None and self.training:
            ei_i, es_i, ev_i = self._dropedge(ei_i, es_i, ev_i, self.dropedge_p)

        h_i, x_i = h0, x
        if ei_i is not None:
            h_skip = h_i
            for i, blk in enumerate(self.blocks_int):
                h_i, x_i = blk(h_i, x_i, ei_i, es_i, ev_i, x_scale=self.x_scale)
                if (i + 1) % 2 == 0 and i > 0:
                    h_i = h_i + 0.2 * h_skip
                    h_skip = h_i

        data_i = data.clone()
        data_i.edge_index = ei_ext
        data_i.edge_attr = (e_s_ext, (e_v_ext if self.use_edge_vec else None))
        if ei_i is not None:
            data_i.internal_edge_index = ei_i
            data_i.internal_edge_attr = (es_i, (ev_i if self.use_edge_vec else None))
        xpool_i, bpool_i = self._hier_pool(data_i, h_i, x_i)
        g_i = self.pool_int(xpool_i, bpool_i)

        # ===== head =====
        g = torch.cat([g_i, g_e], dim=1)  # [B, 192]

        g_main = F.relu(self.fc1(g))
        g_main = F.dropout(g_main, p=self.head_dropout, training=self.training)

        g_res = self.fc_res(g)

        g_combined = g_main + 0.3 * g_res
        g_combined = F.relu(self.fc2(g_combined))
        g_combined = F.dropout(g_combined, p=self.head_dropout * 0.7, training=self.training)

        return self.fc3(g_combined)
