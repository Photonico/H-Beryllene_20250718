"""Independent physical identities and plotting-entry regressions for optics."""

from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from vmatplot import algorithms as algorithms
from vmatplot import linear_optical_properties as optics


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_frequency_units_and_inverse_conversions():
    energy = np.array([0, 0.5, 1, 12.0])
    frequency = algorithms.energy_to_frequency(energy)
    angular = algorithms.energy_to_angular_frequency(energy)
    np.testing.assert_allclose(frequency[2], 1 / 4.135667662e-15)
    np.testing.assert_allclose(angular, 2 * np.pi * frequency, rtol=1e-15)
    np.testing.assert_allclose(algorithms.frequency_to_energy(frequency), energy)
    np.testing.assert_allclose(algorithms.angular_frequency_to_energy(angular), energy)
    assert algorithms.frequency_to_energy(1.0) == algorithms.h_ev


def test_absorption_matches_complex_index_and_vacuum_wavelength():
    energy = np.array([0, 0.5, 2, 6, 12.0])
    epsilon = np.array([2 + 1j, -3 + 0.5j, 0 + 2j, 4 + 3j, 1 + 0j])
    index = np.sqrt(epsilon)
    # Independent optical identity alpha = 4*pi*Im(sqrt(epsilon))/lambda.
    expected = 4 * np.pi * index.imag * energy / (4.135667662e-15 * 2.99792458e17)
    frequency = algorithms.energy_to_frequency(energy)
    actual = optics.comp_absorption_coefficient(frequency, epsilon.real, epsilon.imag)
    np.testing.assert_allclose(actual, expected, atol=1e-15)
    np.testing.assert_allclose(
        optics.current_lop("absorption", frequency, epsilon.real, epsilon.imag), expected
    )
    np.testing.assert_allclose(
        optics.comp_absorption_coefficient_from_angular_frequency(
            algorithms.energy_to_angular_frequency(energy), epsilon.real, epsilon.imag
        ), expected
    )
    legacy = 2 * frequency * index.imag / 2.99792458e17
    nonzero = legacy != 0
    np.testing.assert_allclose(actual[nonzero] / legacy[nonzero], 2 * np.pi)
    assert actual[0] == 0


def test_derived_properties_match_complex_dielectric_identities():
    epsilon = np.array([-3 + 0.5j, 0 + 2j, 4 + 3j, 1 + 0j])
    index = np.sqrt(epsilon)
    n = optics.comp_refractive_index(epsilon.real, epsilon.imag)
    k = optics.comp_extinction_coefficient(epsilon.real, epsilon.imag)
    np.testing.assert_allclose(n, index.real)
    np.testing.assert_allclose(k, index.imag)
    np.testing.assert_allclose(n ** 2 - k ** 2, epsilon.real, atol=1e-14)
    np.testing.assert_allclose(2 * n * k, epsilon.imag, atol=1e-14)
    np.testing.assert_allclose(
        optics.comp_reflectivity(epsilon.real, epsilon.imag),
        np.abs((index - 1) / (index + 1)) ** 2,
    )
    np.testing.assert_allclose(
        optics.comp_energy_loss_spectrum(epsilon.real, epsilon.imag), (-1 / epsilon).imag
    )


@pytest.fixture
def synthetic_system(monkeypatch):
    energy = np.array([0, 1, 2, 4, 8.0])
    raw = {"density_energy_real": energy, "density_energy_imag": energy}
    for position, component in enumerate(["xx", "yy", "zz", "xy", "yx", "yz", "zy", "zx", "xz"]):
        raw[f"density_{component}_real"] = np.array([2, -1, 3, 0, 1.0]) + position / 10
        raw[f"density_{component}_imag"] = np.array([0.5, 2, 1, 3, 0.2])
    data = ["sample", raw, "blue", "solid", 1.5, 1.0, (40, 5)]
    monkeypatch.setattr(optics, "dielectric_systems_list", lambda _: [data])
    return data


@pytest.mark.parametrize("plotter", [
    optics.plot_linear_optical_property,
    optics.plot_linear_optical_property_backup,
    optics.plot_merged_linear_optical_property,
])
@pytest.mark.parametrize("unit,boundary", [("eV", (1, 8)), ("nm", (150, 700))])
@pytest.mark.parametrize("property_name", ["absorption", "refractive", "extinction", "reflectivity", "energy-loss"])
def test_plot_entries_use_normalized_derived_spectra(synthetic_system, plotter, unit, boundary, property_name):
    plotter("test", systems=[["sample", "unused"]], properties=property_name,
            components={"xx": "in-plane", "zz": "out-of-plane"},
            unit=unit, photon_boundary=boundary)
    lines = [line for ax in plt.gcf().axes for line in ax.lines]
    assert len(lines) == 2
    energy = synthetic_system[1]["density_energy_real"]
    x = energy if unit == "eV" else algorithms.energy_to_wavelength(energy)
    selected = np.isfinite(x) & (x >= boundary[0]) & (x <= boundary[1])
    for component, line in zip(["xx", "zz"], lines):
        epsilon = 1 + 8 * (synthetic_system[1][f"density_{component}_real"] - 1)
        epsilon = epsilon + 8j * synthetic_system[1][f"density_{component}_imag"]
        index = np.sqrt(epsilon)
        expected = {
            "absorption": 4 * np.pi * index.imag * energy / (4.135667662e-15 * 2.99792458e17),
            "refractive": index.real,
            "extinction": index.imag,
            "reflectivity": np.abs((index - 1) / (index + 1)) ** 2,
            "energy-loss": (-1 / epsilon).imag,
        }[property_name]
        np.testing.assert_allclose(line.get_xdata(), x[selected])
        np.testing.assert_allclose(line.get_ydata(), expected[selected], atol=1e-14)


@pytest.mark.parametrize("plotter", [optics.plot_linear_optical_property, optics.plot_linear_optical_property_backup])
def test_component_limits_preserve_requested_headroom(synthetic_system, plotter):
    plotter("test", systems=[["sample", "unused"]], properties="energy-loss",
            components=["xx", "zz"], value_boundary={"xx": (-0.1, 2.5), "zz": (-0.1, 9.5)})
    np.testing.assert_allclose(plt.gcf().axes[0].get_ylim(), (-0.1, 2.5))
    np.testing.assert_allclose(plt.gcf().axes[1].get_ylim(), (-0.1, 9.5))


@pytest.mark.parametrize("components", ["xx", ["xx", "yy", "zz", "xy", "yx", "yz", "zy", "zx", "xz"]])
def test_merged_single_and_nine_components(synthetic_system, components):
    optics.plot_merged_linear_optical_property(
        "test", systems=[["sample", "unused"]], properties="absorption", components=components
    )
    assert len(plt.gca().lines) == (1 if isinstance(components, str) else 9)


@pytest.mark.parametrize("material,thickness,expected", [
    ("a-Beryllene", 3.96, [0.135476, 0.175010]),
    ("b-Beryllene", 5.85508656790224, [0.159298, 0.241453]),
])
def test_archived_pristine_absorption_from_density_density_hdf5(material, thickness, expected):
    path = Path(__file__).resolve().parents[1] / "6.0_dielectric_selection" / material / "vaspout.h5"
    with h5py.File(path, "r") as handle:
        energy = handle["results/linear_response/energies_dielectric_function"][:]
        tensor = handle["results/linear_response/density_density_dielectric_function"][:]
    selected = (energy >= 0) & (energy <= 12)
    maxima = []
    for component in [0, 2]:
        epsilon = 1 + (40 / thickness) * (tensor[component, component, :, 0] - 1)
        epsilon = epsilon + 1j * (40 / thickness) * tensor[component, component, :, 1]
        alpha = optics.comp_absorption_coefficient(
            algorithms.energy_to_frequency(energy), epsilon.real, epsilon.imag
        )
        independent = 4 * np.pi * np.sqrt(epsilon).imag * energy / (4.135667662e-15 * 2.99792458e17)
        np.testing.assert_allclose(alpha, independent, atol=1e-12)
        maxima.append(alpha[selected].max())
    np.testing.assert_allclose(maxima, expected, atol=0.0000005, rtol=0)
