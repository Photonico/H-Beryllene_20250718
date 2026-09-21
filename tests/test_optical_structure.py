"""Regression tests for optical thickness, units, and periodic geometry."""

from pathlib import Path

import numpy as np
import pytest

from vmatplot.structure import slab_optical_factor, slab_optical_geometry


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("name, expected, surfaces", [
    ("a-Beryllene", 3.96, ("Be", "Be")),
    ("b-Beryllene", 5.85508656790224, ("Be", "Be")),
    ("c3-Beryllene", 6.966489780213428, ("Be", "Be")),
    ("a-Beryllene_hh", 3.144512351252184, ("H", "H")),
    ("b-Beryllene_b", 5.523544103908468, ("Be", "H")),
    ("b-Beryllene_bb", 5.192831416812704, ("H", "H")),
])
def test_selected_structure_thickness(name, expected, surfaces):
    geometry = slab_optical_geometry(ROOT / "6.0_dielectric_selection" / name)
    assert geometry.coordinate_mode == "Direct"
    assert geometry.cell_height == pytest.approx(40)
    assert geometry.effective_thickness == pytest.approx(expected, abs=1e-12)
    assert (geometry.lower_element, geometry.upper_element) == surfaces
    assert not geometry.crosses_boundary


def write_structure(path, lattice, positions, mode="Direct", scale="1", selective=False):
    lines = ["Synthetic isolated slab", scale]
    lines += [" ".join(map(str, row)) for row in lattice]
    lines += ["Be H", "1 1"]
    if selective:
        lines += ["Selective dynamics"]
    lines += [mode]
    lines += [" ".join(map(str, row)) + (" T T F" if selective else "") for row in positions]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_periodic_translation_does_not_change_thickness(tmp_path):
    lattice = np.diag([3, 3, 40])
    centred = write_structure(tmp_path / "centred", lattice, [[0, 0, .1], [0, 0, .2]])
    wrapped = write_structure(tmp_path / "wrapped", lattice, [[0, 0, .95], [0, 0, .05]])
    assert slab_optical_factor(centred) == pytest.approx((40, 6.68))
    assert slab_optical_factor(wrapped) == pytest.approx(slab_optical_factor(centred))
    assert slab_optical_geometry(wrapped).crosses_boundary


def test_scaled_cartesian_matches_direct_with_selective_dynamics(tmp_path):
    lattice = np.diag([3, 3, 20])
    direct = write_structure(tmp_path / "direct", lattice, [[0, 0, .1], [0, 0, .2]], scale="2")
    cartesian = write_structure(tmp_path / "cartesian", lattice, [[0, 0, 2], [0, 0, 4]],
                                mode="Cartesian", scale="2", selective=True)
    assert slab_optical_factor(cartesian) == pytest.approx(slab_optical_factor(direct))
    assert slab_optical_geometry(cartesian).coordinate_mode == "Cartesian"


def test_cell_height_is_normal_projection_and_rotation_invariant(tmp_path):
    # The third vector length exceeds 40 A, but the normal repeat is exactly 40 A.
    lattice = np.array([[3., 0, 0], [0, 3, 0], [10, 0, 40]])
    fractions = np.array([[0, 0, .1], [0, 0, .2]])
    angle = np.pi / 3
    rotation = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    path = write_structure(tmp_path / "rotated", lattice @ rotation, fractions)
    assert slab_optical_factor(path) == pytest.approx((40, 6.68))


def test_negative_scale_is_requested_volume(tmp_path):
    lattice = np.diag([1, 1, 10])
    path = write_structure(tmp_path / "negative", lattice, [[0, 0, .1], [0, 0, .2]], scale="-640")
    assert slab_optical_factor(path) == pytest.approx((40, 6.68))


def test_surface_radius_convention_is_not_atomic_sphere_envelope():
    # The retained H/H endpoint convention is intentionally thinner than 2*r_Be.
    geometry = slab_optical_geometry(ROOT / "6.0_dielectric_selection/a-Beryllene_hh")
    assert geometry.effective_thickness < 2 * 1.98
    assert geometry.lower_radius == geometry.upper_radius == .70


def test_bulk_sized_effective_thickness_is_rejected(tmp_path):
    path = write_structure(tmp_path / "bulk", np.diag([3, 3, 2]), [[0, 0, 0], [0, 0, .5]])
    with pytest.raises(ValueError, match="below the cell height"):
        slab_optical_factor(path)
