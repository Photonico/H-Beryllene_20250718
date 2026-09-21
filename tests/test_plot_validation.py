"""Regression checks for visible optical peaks and x-window intersections."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from vmatplot.plot_validation import assert_curve_visibility, ensure_curve_visibility


@pytest.fixture
def axes():
    fig, ax = plt.subplots()
    yield fig, ax
    plt.close(fig)


@pytest.mark.parametrize("peak,limits", [(8.878109, (-0.1, 2.1)), (-8.0, (-2.1, 0.1))])
def test_visible_clipped_peak_rejected_then_included(axes, peak, limits):
    fig, ax = axes
    ax.plot([0, 1, 2], [0, peak, 0], label="loss peak")
    ax.set_xlim(0, 2)
    ax.set_ylim(*limits)
    with pytest.raises(ValueError, match="loss peak"):
        assert_curve_visibility(fig)
    assert ensure_curve_visibility(fig) is fig
    low, high = ax.get_ylim()
    assert low <= limits[0] and high >= limits[1]
    assert low < peak < high
    assert assert_curve_visibility(fig) == 1


def test_outside_window_values_do_not_expand_axis(axes):
    fig, ax = axes
    ax.plot([-2, -1, 0, 1, 2, 3], [1000, -1000, 1, 2, 1000, -1000], label="spectrum")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 3)
    assert assert_curve_visibility(fig) == 1
    ensure_curve_visibility(fig)
    np.testing.assert_allclose(ax.get_ylim(), (0, 3))


def test_crossing_segment_checked_without_any_samples_in_window(axes):
    fig, ax = axes
    # This line is y = 10*x. Its visible endpoints are (0, 0) and (1, 10),
    # although neither stored sample lies inside the visible x window.
    ax.plot([-1, 2], [-10, 20], label="crossing spectrum")
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, 9)
    with pytest.raises(ValueError, match="crossing spectrum"):
        assert_curve_visibility(fig)
    ensure_curve_visibility(fig)
    low, high = ax.get_ylim()
    assert -10 < low < 0
    assert 10 < high < 20
    assert assert_curve_visibility(fig) == 1


def test_visible_spectrum_background_follows_expanded_y_limits(axes):
    fig, ax = axes
    background = ax.imshow(np.zeros((2, 2)), extent=(0.5, 1.5, -0.1, 2.1), aspect="auto")
    ax.plot([0, 1, 2], [0, 8.878109, 0], label="loss peak")
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.1, 2.1)
    ensure_curve_visibility(fig)
    np.testing.assert_allclose(background.get_extent()[:2], (0.5, 1.5))
    np.testing.assert_allclose(background.get_extent()[2:], ax.get_ylim())
    np.testing.assert_allclose(ax.get_xlim(), (0, 2))
