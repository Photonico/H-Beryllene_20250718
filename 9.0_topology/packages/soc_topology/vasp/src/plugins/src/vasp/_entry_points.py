import importlib
import warnings
from importlib.metadata import entry_points

try:
    vasp_plugin = importlib.import_module("vasp_plugin")
except ModuleNotFoundError:
    vasp_plugin = None


def load(function_name):
    functions = [
        entry_point.load()
        for entry_point in entry_points(name=function_name, group="vasp")
    ]
    if hasattr(vasp_plugin, function_name):
        functions.append(getattr(vasp_plugin, function_name))
    if not functions:
        message = """\
No plugin for {} found. Please make sure that you have either installed Python packages \
that provide entry points for VASP, installed a package called vasp_plugin, or have a \
file called vasp_plugin.py in the directory in which you have launched your VASP \
calculation."""
        warnings.warn(message.format(function_name))
    return functions
