"""Structure-derived factors for the project's effective slab response.

The retained convention is the span of the outermost atomic *centres* plus
the radii of the atoms at those two surfaces (Be: 1.98 A; H: 0.70 A).
This is a chosen optical thickness, not a union of atomic spheres or a
unique physical thickness. Bulk calculations must keep the factor (1, 1).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


OPTICAL_RADII_ANGSTROM = {"Be": 1.98, "H": 0.70}


@dataclass(frozen=True)
class SlabOpticalGeometry:
    """Lengths are in angstrom; fractions refer to the surface-normal period."""

    structure_file: str
    coordinate_mode: str
    cell_height: float
    fractional_span: float
    atomic_span: float
    lower_element: str
    upper_element: str
    lower_radius: float
    upper_radius: float
    effective_thickness: float
    crosses_boundary: bool

    @property
    def factor(self):
        """Return (full supercell height, effective slab thickness)."""
        return self.cell_height, self.effective_thickness


def _read_vasp_structure(path):
    """Read VASP 5/6 positions, including scale and selective dynamics."""
    path = Path(path)
    if path.is_dir():
        # Use the relaxed geometry; POSCAR is only a fallback if unavailable.
        path = next((path / name for name in ("CONTCAR", "POSCAR")
                     if (path / name).is_file()), path / "CONTCAR")
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        raise ValueError(f"Incomplete VASP structure: {path}")
    scale_values = np.asarray(lines[1].split(), dtype=float)
    raw_lattice = np.asarray([line.split()[:3] for line in lines[2:5]], dtype=float)
    if scale_values.size == 1:
        scale = float(scale_values[0])
        if scale < 0:
            # A negative scalar specifies the requested cell volume.
            scale = (abs(scale) / abs(np.linalg.det(raw_lattice))) ** (1 / 3)
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid VASP scaling factor: {path}")
    elif scale_values.size == 3 and np.all(scale_values > 0):
        scale = scale_values
    else:
        raise ValueError(f"Expected one nonzero or three positive scale factors: {path}")
    lattice = raw_lattice * scale
    elements = lines[5].split()
    try:
        counts = [int(value) for value in lines[6].split()]
    except ValueError as exc:
        raise ValueError(f"Element labels are required (VASP 5/6 format): {path}") from exc
    if len(elements) != len(counts) or not counts or any(n <= 0 for n in counts):
        raise ValueError(f"Invalid element counts: {path}")
    line_index = 7
    if lines[line_index].strip().lower().startswith("s"):
        line_index += 1
    mode = lines[line_index].strip().lower()
    if mode.startswith("d"):
        coordinate_mode = "Direct"
    elif mode.startswith(("c", "k")):
        coordinate_mode = "Cartesian"
    else:
        raise ValueError(f"Unrecognised VASP coordinate mode: {path}")
    natoms = sum(counts)
    positions = np.asarray([
        line.split()[:3] for line in lines[line_index + 1:line_index + 1 + natoms]
    ], dtype=float)
    if positions.shape != (natoms, 3) or not np.all(np.isfinite(positions)):
        raise ValueError(f"Incomplete or nonfinite atomic positions: {path}")
    if not np.all(np.isfinite(lattice)) or abs(np.linalg.det(lattice)) < 1e-12:
        raise ValueError(f"Nonfinite or singular lattice: {path}")
    cartesian = positions @ lattice if coordinate_mode == "Direct" else positions * scale
    species = np.repeat(elements, counts)
    return path, coordinate_mode, lattice, cartesian, species


def slab_optical_geometry(path, radii=None):
    """Read an isolated slab and apply the original surface-radius convention.

    The slab plane is spanned by the first two lattice vectors. The full cell
    height is the projection of the third vector onto its normal, rather than
    its length or a user-supplied vacuum thickness. The largest cyclic gap in
    projected atomic centres identifies the vacuum; cutting there handles a
    single slab crossing the periodic boundary. This routine is intentionally
    not a bulk/slab classifier and assumes one isolated slab per cell.
    """
    radii = OPTICAL_RADII_ANGSTROM if radii is None else radii
    path, mode, lattice, positions, species = _read_vasp_structure(path)
    normal = np.cross(lattice[0], lattice[1])
    normal /= np.linalg.norm(normal)
    signed_height = float(np.dot(lattice[2], normal))
    if signed_height < 0:
        normal = -normal
    cell_height = abs(signed_height)
    fractions = np.mod(positions @ normal / cell_height, 1.0)
    sorted_fractions = np.sort(fractions)
    gaps = np.diff(np.concatenate((sorted_fractions, sorted_fractions[:1] + 1)))
    vacuum_index = int(np.argmax(gaps))
    start = sorted_fractions[(vacuum_index + 1) % len(fractions)]
    unwrapped = np.mod(fractions - start, 1.0)
    # Floating point subtraction can put an atom at the cut just below 1.
    unwrapped[np.isclose(unwrapped, 1.0, rtol=0, atol=1e-12)] = 0.0
    fractional_span = float(np.max(unwrapped))
    lower_species = species[np.isclose(unwrapped, 0.0, rtol=0, atol=1e-10)]
    upper_species = species[np.isclose(unwrapped, fractional_span, rtol=0, atol=1e-10)]
    try:
        # For tied surface positions, use the largest specified endpoint radius.
        lower_element = max(lower_species, key=lambda element: radii[element])
        upper_element = max(upper_species, key=lambda element: radii[element])
        lower_radius = float(radii[lower_element])
        upper_radius = float(radii[upper_element])
    except KeyError as exc:
        raise ValueError(f"Missing optical radius for {exc.args[0]} in {path}") from exc
    if not all(np.isfinite(r) and r > 0 for r in (lower_radius, upper_radius)):
        raise ValueError("Optical radii must be finite and positive")
    atomic_span = fractional_span * cell_height
    thickness = atomic_span + lower_radius + upper_radius
    if thickness >= cell_height:
        raise ValueError(f"Effective slab thickness must be below the cell height: {path}")
    return SlabOpticalGeometry(
        str(path), mode, cell_height, fractional_span, atomic_span,
        str(lower_element), str(upper_element), lower_radius, upper_radius,
        thickness, bool(np.ptp(fractions) > fractional_span + 1e-10),
    )


def slab_optical_factor(path, radii=None):
    """Return the dielectric normalization factor from this structure."""
    return slab_optical_geometry(path, radii=radii).factor
