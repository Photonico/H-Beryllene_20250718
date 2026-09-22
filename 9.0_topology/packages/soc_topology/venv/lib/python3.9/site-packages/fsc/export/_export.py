#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    08.04.2016 11:43:47 CEST
# File:    _export.py
"""
Defines a decorator / function for adding object to the current
module's __all__. Optionally, a check can be enabled which ensures
that exported objects have a docstring.
"""

from typing import Optional, TypeVar

import importlib

__all__ = ['test_doc', 'export']

enable_doc_test: bool = False


def test_doc(enable: bool = True) -> None:
    """Test whether exported objects have a non-zero docstring. This must be called before the module is imported.

    :param enable:  Determines whether the test should be enabled or disabled.
    """
    global enable_doc_test  # pylint: disable=global-statement,invalid-name
    enable_doc_test = enable


T = TypeVar('T')  # pylint: disable=invalid-name


def export(obj: T, name: Optional[str] = None) -> T:
    """Adds the decorated object to its module's ``__all__``.

    :param name:    Optional parameter to overwrite the object's name that is written to ``__all__``.
    :type name: str
    """
    obj_name: str
    if name is None:
        obj_name = obj.__name__  # type: ignore
    else:
        obj_name = name

    mod = importlib.import_module(obj.__module__)
    # add to __all__
    try:
        mod.__all__.append(obj_name)  # type: ignore
    except AttributeError:
        mod.__all__ = [obj_name]  # type: ignore

    # test if the docstring is nonzero
    if enable_doc_test:
        if not obj.__doc__:
            raise AssertionError("'{}' does not have a non-zero docstring.".format(obj_name))

    return obj
