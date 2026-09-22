from unittest.mock import MagicMock, patch

import pytest

import vasp._entry_points as entry_points


def test_nonexisting_plugin(function_name):
    with pytest.warns(UserWarning):
        assert entry_points.load(function_name) == []


@patch("vasp._entry_points.vasp_plugin")
def test_vasp_plugin(mock_plugin, function_name):
    function = getattr(mock_plugin, function_name)
    assert entry_points.load(function_name) == [function]
    delattr(mock_plugin, function_name)
    with pytest.warns(UserWarning):
        assert entry_points.load(function_name) == []


@patch("vasp._entry_points.entry_points", return_value=(MagicMock(), MagicMock()))
def test_entry_points(mock_entry_points, function_name):
    functions = mock_entry_points.return_value
    expected = [function.load.return_value for function in functions]
    assert entry_points.load(function_name) == expected
    for function in functions:
        function.load.assert_called_once_with()
        function.reset_mock()
