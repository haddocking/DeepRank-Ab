"""Unit tests for pure math helpers in src/tools/edge_orientation.py (CPU tensors only)."""

import torch

from src.tools.edge_orientation import (
    compute_edge_orientation,
    construct_3d_basis,
    normalize_vector,
    project_v2v,
)


def test_normalize_vector_unit_length():
    v = torch.tensor([[3.0, 4.0, 0.0]])
    normed = normalize_vector(v, dim=-1)
    assert torch.allclose(normed, torch.tensor([[0.6, 0.8, 0.0]]), atol=1e-4)
    assert torch.allclose(normed.norm(dim=-1), torch.tensor([1.0]), atol=1e-4)


def test_normalize_vector_handles_zero_vector():
    v = torch.zeros((1, 3))
    normed = normalize_vector(v, dim=-1)
    assert torch.allclose(normed, torch.zeros((1, 3)))


def test_project_v2v_onto_unit_axis():
    v = torch.tensor([[1.0, 2.0, 3.0]])
    e = torch.tensor([[1.0, 0.0, 0.0]])  # unit vector along x
    proj = project_v2v(v, e, dim=-1)
    assert torch.allclose(proj, torch.tensor([[1.0, 0.0, 0.0]]))


def test_project_v2v_orthogonal_is_zero():
    v = torch.tensor([[0.0, 1.0, 0.0]])
    e = torch.tensor([[1.0, 0.0, 0.0]])
    proj = project_v2v(v, e, dim=-1)
    assert torch.allclose(proj, torch.zeros((1, 3)))


def test_construct_3d_basis_orthonormal():
    center = torch.tensor([[0.0, 0.0, 0.0]])
    p1 = torch.tensor([[1.0, 0.0, 0.0]])
    p2 = torch.tensor([[0.0, 1.0, 0.0]])

    basis = construct_3d_basis(center, p1, p2)  # (1, 3, 3), columns [e1, e2, e3]
    assert basis.shape == (1, 3, 3)

    e1, e2, e3 = basis[0, :, 0], basis[0, :, 1], basis[0, :, 2]
    assert torch.allclose(e1, torch.tensor([1.0, 0.0, 0.0]), atol=1e-4)
    assert torch.allclose(e2, torch.tensor([0.0, 1.0, 0.0]), atol=1e-4)
    assert torch.allclose(e3, torch.tensor([0.0, 0.0, 1.0]), atol=1e-4)

    # basis vectors should be mutually orthogonal and unit length
    for a in (e1, e2, e3):
        assert torch.allclose(a.norm(), torch.tensor(1.0), atol=1e-4)
    assert torch.allclose(torch.dot(e1, e2), torch.tensor(0.0), atol=1e-4)
    assert torch.allclose(torch.dot(e1, e3), torch.tensor(0.0), atol=1e-4)


def test_compute_edge_orientation_single_edge():
    # Two nodes sharing an identical local frame (e1=x, e2=y, e3=z), offset along x.
    pos_CA = torch.tensor([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    pos_C = torch.tensor([[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    pos_N = torch.tensor([[0.0, 1.0, 0.0], [2.0, 1.0, 0.0]])
    edge_index = torch.tensor([[0], [1]])  # src=0, dst=1

    orientation = compute_edge_orientation(pos_CA, pos_C, pos_N, edge_index)

    assert orientation.shape == (1, 3)
    # direction src->dst is -x in world space, and src's local frame is identity,
    # so expressed direction should be (-1, 0, 0).
    assert torch.allclose(orientation, torch.tensor([[-1.0, 0.0, 0.0]]), atol=1e-4)
