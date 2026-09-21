"""Check the actual plotted samples before exporting scientific figures."""

import matplotlib.pyplot as plt
import numpy as np


def visible_curve_values(ax):
    """Yield finite samples and interpolated boundary values inside the x window."""
    left, right = sorted(ax.get_xlim())
    for line in ax.lines:
        if not line.get_visible():
            continue
        x, y = np.asarray(line.get_xdata(), float), np.asarray(line.get_ydata(), float)
        finite = np.isfinite(x) & np.isfinite(y)
        values = y[finite & (x >= left) & (x <= right)]
        boundary_values = []
        for boundary in (left, right):
            crossing = finite[:-1] & finite[1:] & ((x[:-1] - boundary) * (x[1:] - boundary) < 0)
            indices = np.flatnonzero(crossing)
            boundary_values.extend(y[indices] + (boundary - x[indices]) *
                                   (y[indices + 1] - y[indices]) / (x[indices + 1] - x[indices]))
        if boundary_values:
            values = np.concatenate([values, boundary_values])
        if values.size:
            yield line, values


def ensure_curve_visibility(fig=None, margin=0.05):
    """Expand, never shrink, each y range to include every visible curve sample."""
    fig = plt.gcf() if fig is None else fig
    if margin < 0:
        raise ValueError("margin must be non-negative")
    for ax in fig.axes:
        arrays = [values for _, values in visible_curve_values(ax)]
        if not arrays:
            continue
        low, high = ax.get_ylim()
        if low >= high:
            raise ValueError("Optical figure validation requires increasing y limits")
        values = np.concatenate(arrays)
        minimum, maximum = float(values.min()), float(values.max())
        padding = margin * max(maximum - minimum, abs(maximum), abs(minimum), 1e-12)
        updated = (min(low, minimum - padding), max(high, maximum + padding))
        ax.set_ylim(*updated)
        # Preserve the visible-spectrum background when the range expands.
        for background in ax.images:
            x0, x1, _, _ = background.get_extent()
            background.set_extent((x0, x1, *updated))
        ax.set_ylim(*updated)
    assert_curve_visibility(fig)
    return fig


def assert_curve_visibility(fig=None):
    """Raise if any finite plotted sample in the x window is clipped in y."""
    fig = plt.gcf() if fig is None else fig
    checked = 0
    for index, ax in enumerate(fig.axes):
        low, high = sorted(ax.get_ylim())
        tolerance = 1e-10 * max(abs(low), abs(high), 1)
        for line, values in visible_curve_values(ax):
            if values.min() < low - tolerance or values.max() > high + tolerance:
                raise ValueError(f"Axis {index}, {line.get_label()}: data [{values.min()}, "
                                 f"{values.max()}] outside y limits [{low}, {high}]")
            checked += 1
    return checked
