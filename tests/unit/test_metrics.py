"""Unit tests for pure numpy/sklearn math in src/Metrics.py.

Note: `_safe_div` has a real bug for array inputs — `np.divide(num, den,
where=(den != 0))` is called without an `out=` array, so positions where the
denominator is zero contain uninitialized memory (observed: values like
1.14e-313), not 0 as the "zero-denominator protection" docstring implies.
This only affects the array path (scalar calls happen to come back as 0.0
because of how numpy allocates 0-d arrays) — flagged for the user rather than
asserted on here, since asserting on undefined memory would be meaningless.
"""

import numpy as np
import pytest
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score

from src.Metrics import (
    Metrics,
    _finite_pairmask,
    _safe_div,
    _to_np,
    get_binary,
    get_comparison,
)


# ---------------------------------------------------------------------------
# _to_np
# ---------------------------------------------------------------------------

def test_to_np_converts_list_to_float_array():
    out = _to_np([1, 2, 3])
    assert out.dtype == np.float64
    np.testing.assert_array_equal(out, [1.0, 2.0, 3.0])


# ---------------------------------------------------------------------------
# _finite_pairmask
# ---------------------------------------------------------------------------

def test_finite_pairmask_filters_nans():
    y, pred, mask = _finite_pairmask([1, float("nan"), 3], [1, 2, float("nan")])
    np.testing.assert_array_equal(y, [1.0])
    np.testing.assert_array_equal(pred, [1.0])
    np.testing.assert_array_equal(mask, [True, False, False])


def test_finite_pairmask_all_finite():
    y, pred, mask = _finite_pairmask([1, 2, 3], [4, 5, 6])
    assert mask.all()
    np.testing.assert_array_equal(y, [1.0, 2.0, 3.0])


# ---------------------------------------------------------------------------
# _safe_div
# ---------------------------------------------------------------------------

def test_safe_div_scalar_nonzero():
    assert _safe_div(6, 3) == 2.0


def test_safe_div_scalar_zero_denominator_does_not_raise():
    # Documented behavior: division by zero is guarded (no ZeroDivisionError/warning).
    result = _safe_div(5, 0)
    assert result == 0.0


def test_safe_div_array_nonzero_denominators():
    out = _safe_div(np.array([2.0, 4.0]), np.array([2.0, 2.0]))
    np.testing.assert_allclose(out, [1.0, 2.0])


# ---------------------------------------------------------------------------
# get_binary
# ---------------------------------------------------------------------------

def test_get_binary_dockq_is_greater_than_threshold():
    # dockq is in the "inverse" set: higher value -> positive class
    assert get_binary([0.1, 0.5, 0.9], threshold=0.23, target="dockq") == [0, 1, 1]


def test_get_binary_irmsd_is_less_than_threshold():
    # irmsd is not in the inverse set: lower value -> positive class
    assert get_binary([0.1, 0.5, 0.9], threshold=0.5, target="irmsd") == [1, 0, 0]


# ---------------------------------------------------------------------------
# get_comparison
# ---------------------------------------------------------------------------

def test_get_comparison_binary_counts():
    prediction = [1, 0, 1, 1]
    ground_truth = [1, 0, 0, 1]
    fp, fn, tp, tn = get_comparison(prediction, ground_truth, binary=True)
    assert (fp, fn, tp, tn) == (1, 0, 2, 1)


def test_get_comparison_handles_missing_class_gracefully():
    # only class 0 present -> should not raise, class 1 counts are 0
    prediction = [0, 0, 0]
    ground_truth = [0, 0, 0]
    fp, fn, tp, tn = get_comparison(prediction, ground_truth, binary=True)
    assert (fp, fn, tp) == (0, 0, 0)


# ---------------------------------------------------------------------------
# Metrics — regression branch
# ---------------------------------------------------------------------------

def test_metrics_regression_computes_sklearn_scores():
    prediction = [0.1, 0.6, 0.9, 0.4]
    y = [0.05, 0.55, 0.95, 0.35]
    m = Metrics(prediction, y, target="dockq", threshold=0.5)

    assert m.mean_absolute_error == pytest.approx(mean_absolute_error(y, prediction))
    assert m.r2_score == pytest.approx(r2_score(y, prediction))
    assert m.root_mean_squared_error == pytest.approx(m.mean_squared_error ** 0.5)


def test_metrics_regression_binary_confusion_counts():
    # threshold=0.5, dockq inverse convention (value > threshold -> positive)
    prediction = [0.9, 0.1, 0.8, 0.2]
    y = [0.9, 0.05, 0.05, 0.9]
    m = Metrics(prediction, y, target="dockq", threshold=0.5)
    # predicted positives: idx 0, 2 ; true positives among those: idx 0 only
    assert m.precision == pytest.approx(0.5)


def test_metrics_empty_input_leaves_all_none():
    m = Metrics([], [], target="dockq", threshold=0.5)
    assert m.sensitivity is None
    assert m.accuracy is None
    assert m.r2_score is None


def test_metrics_filters_non_finite_pairs():
    prediction = [0.1, float("nan"), 0.9]
    y = [0.05, 0.5, 0.95]
    m = Metrics(prediction, y, target="dockq", threshold=0.5)
    assert m.prediction.size == 2
    assert m.y.size == 2


# ---------------------------------------------------------------------------
# Metrics — ranking utilities
# ---------------------------------------------------------------------------

@pytest.fixture
def ranking_metrics():
    prediction = [0.9, 0.1, 0.8, 0.2]
    y = [0.9, 0.05, 0.05, 0.9]
    return Metrics(prediction, y, target="dockq", threshold=0.5)


def test_format_score_reverses_order_for_inverse_targets(ranking_metrics):
    idx, ground_truth_bool = ranking_metrics.format_score()
    # dockq is "inverse": highest prediction ranked first
    np.testing.assert_array_equal(idx, [0, 2, 3, 1])
    np.testing.assert_array_equal(ground_truth_bool, [1, 0, 0, 1])


def test_hitrate_cumulative_sum(ranking_metrics):
    hits = ranking_metrics.hitrate()
    np.testing.assert_array_equal(hits, [1, 1, 2, 2])


def test_hitrate_at_k(ranking_metrics):
    assert ranking_metrics.hitrate_at_k(k=2) == pytest.approx(0.5)


def test_success_at_k(ranking_metrics):
    assert ranking_metrics.success_at_k(k=2) == 1.0


def test_success_at_k_no_hits_in_topk():
    prediction = [0.9, 0.1, 0.8, 0.2]
    y = [0.05, 0.9, 0.05, 0.9]  # top-2 ranked predictions (idx 0, 2) are both negatives
    m = Metrics(prediction, y, target="dockq", threshold=0.5)
    assert m.success_at_k(k=2) == 0.0


def test_auc_matches_sklearn(ranking_metrics):
    y_binary = get_binary([0.9, 0.05, 0.05, 0.9], threshold=0.5, target="dockq")
    expected = roc_auc_score(y_binary, [0.9, 0.1, 0.8, 0.2])
    assert ranking_metrics.auc() == pytest.approx(expected)


def test_auc_none_when_single_class():
    prediction = [0.1, 0.2, 0.3]
    y = [0.05, 0.1, 0.15]  # all below threshold -> single class after binarization
    m = Metrics(prediction, y, target="dockq", threshold=0.5)
    assert m.auc() is None
