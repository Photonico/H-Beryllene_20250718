import functools
from dataclasses import dataclass, field
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

import vasp
import vasp._machine_learning
from vasp.typing import IndexArray


@dataclass(frozen=True)
class MockConstants:
    ion_types: IndexArray = field(default_factory=lambda: np.arange(1, 4))


class MockCalculator:
    def __init__(self, shape):
        self._shape = shape

    def get_potential_energy(self, *args, **kwargs):
        return 1.0

    def get_forces(self, *args, **kwargs):
        return np.ones(self._shape)

    def get_stress(self, *args, **kwargs):
        return np.ones((3, 3))


@pytest.fixture
def mock_interface(function_name):
    interface = getattr(vasp, f"interface_{function_name}")
    return functools.partial(interface, MockConstants(), MagicMock())


@pytest.fixture
def mock_plugin():
    with patch("vasp._entry_points.vasp_plugin") as mock:
        yield mock


@pytest.fixture
def mock_entry_points():
    return_value = (MagicMock(), MagicMock())
    with patch("vasp._entry_points.entry_points", return_value=return_value) as mock:
        yield mock


def test_interface(mock_interface):
    with pytest.warns(UserWarning):
        assert mock_interface() == 0


def test_vasp_plugin(mock_plugin, mock_interface, function_name):
    function = getattr(mock_plugin, function_name)
    assert mock_interface() == 0
    function.assert_called_once_with(*mock_interface.args)


def test_vasp_plugin(mock_plugin, mock_interface, function_name):
    function = getattr(mock_plugin, function_name)
    delattr(mock_plugin, function_name)
    with pytest.warns(UserWarning):
        assert mock_interface() == 0
    function.assert_not_called()


@patch("vasp._entry_points.entry_points", return_value=(MagicMock(), MagicMock()))
def test_entry_points(mock_entry_points, mock_interface):
    plugins = mock_entry_points.return_value
    functions = [plugin.load.return_value for plugin in plugins]
    assert mock_interface() == 0
    for function in functions:
        function.assert_called_once_with(*mock_interface.args)
        function.reset_mock()


def test_interface_force_and_stress(mock_plugin):
    constants = vasp.ConstantsForceAndStress(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        forces=np.random.random((3, 5)),
        stress=np.random.random((3, 3)),
    )
    assert vasp.interface_force_and_stress(constants, None) == vasp.SUCCESS
    assert all(constants.ion_types == [0, 1, 2, 2, 2])  # output: Python style indexing
    array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,POMASS,forces,stress,charge_density"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_machine_learning(mock_plugin):
    vasp._machine_learning.get_calculators.cache_clear()
    constants = mock_force_and_stress_constants()
    additions = vasp.AdditionsForceAndStress(
        total_energy=2.0,
        forces=np.full((3, 5), 3.0),
        stress=np.full((3, 3), 4.0),
    )
    mock_plugin.machine_learning.return_value = MockCalculator(constants.forces.shape)
    assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
    assert all(constants.ion_types == [0, 1, 2, 2, 2])  # output: Python style indexing
    array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,POMASS,forces,stress,charge_density"""
    check_arrays_not_writeable(constants, array_names)
    assert np.allclose(additions.total_energy, 3.0)
    assert np.allclose(additions.forces, 4.0)
    assert np.allclose(additions.stress, 5.0)
    # check calling the interface for a second time does not create a new calculator
    constants = mock_force_and_stress_constants()
    assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
    mock_plugin.machine_learning.assert_called_once()
    check_arrays_not_writeable(constants, array_names)
    assert np.allclose(additions.total_energy, 4.0)
    assert np.allclose(additions.forces, 5.0)
    assert np.allclose(additions.stress, 6.0)


def test_entry_points_machine_learning(mock_entry_points):
    vasp._machine_learning.get_calculators.cache_clear()
    constants = mock_force_and_stress_constants()
    additions = vasp.AdditionsForceAndStress(
        total_energy=2.0,
        forces=np.full((3, 5), 3.0),
        stress=np.full((3, 3), 4.0),
    )
    shape = constants.forces.shape
    first_plugin, second_plugin = mock_entry_points.return_value
    first_interface = first_plugin.load.return_value
    first_interface.return_value = MockCalculator(shape)
    second_interface = second_plugin.load.return_value
    second_interface.return_value = MockCalculator(shape)
    assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
    first_interface.assert_called_once()
    second_interface.assert_called_once()
    assert np.allclose(additions.total_energy, 4.0)
    assert np.allclose(additions.forces, 5.0)
    assert np.allclose(additions.stress, 6.0)


def test_interface_local_potential(mock_plugin):
    constants = vasp.ConstantsLocalPotential(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
    )
    assert vasp.interface_local_potential(constants, None) == vasp.SUCCESS
    array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,charge_density"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_structure(mock_plugin):
    shape_grid = [60, 100, 80]
    constants = vasp.ConstantsStructure(
        **mock_SrTiO3_structure(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        total_energy=np.random.random(),
        forces=np.random.random((3, 5)),
        stress=np.random.random((3, 3)),
        shape_grid=shape_grid,
        charge_density=np.random.random(shape_grid),
    )
    assert vasp.interface_structure(constants, None) == vasp.SUCCESS
    array_names = """ion_types,atomic_numbers,lattice_vectors,positions,POMASS,\
forces,stress"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_occupancies(mock_plugin):
    constants = vasp.ConstantsOccupancies(
        NELECT=27.0, EFERMI=2.5, NUPDOWN=0.5, ISMEAR=1, SIGMA=0.1, EMIN=-5.2, EMAX=6.3
    )
    assert vasp.interface_occupancies(constants, None) == vasp.SUCCESS


def test_error_in_plugin(function_name, mock_plugin):
    def broken_interface(*args):
        raise NotImplementedError("There is no implementation here.")

    setattr(mock_plugin, function_name, broken_interface)
    func = getattr(vasp, f"interface_{function_name}")
    assert func(MockConstants(), None) == vasp.ERROR_IN_PLUGIN


def mock_force_and_stress_constants():
    return vasp.ConstantsForceAndStress(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        forces=np.zeros((3, 5)),
        stress=np.zeros((3, 3)),
    )


def mock_SrTiO3_structure():
    return {
        "number_ions": 5,
        "number_ion_types": 3,
        "ion_types": np.array([1, 2, 3, 3, 3]),  # input: Fortran style indexing
        "atomic_numbers": np.array([38, 22, 8, 8, 8]),
        "lattice_vectors": 4.0 * np.eye(3),
        "positions": np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ),
    }


def mock_electron_data():
    shape = (10, 12, 14)
    return {
        "NELECT": 40.0,
        "ENCUT": 250.0,
        "ZVAL": np.array([10.0, 12.0, 6.0]),
        "shape_grid": np.array(shape),
        "charge_density": np.random.random(shape),
    }


def check_arrays_not_writeable(constants, array_names):
    for array_name in array_names.split(","):
        array = getattr(constants, array_name)
        assert not array.flags.writeable
