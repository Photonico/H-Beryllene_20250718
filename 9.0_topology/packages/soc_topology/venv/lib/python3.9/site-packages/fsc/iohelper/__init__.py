#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  C. Frescolino, D. Gresch
# File:    __init__.py

"""
This is a helper for saving and loading data from files, with a given encoding and decoding function. It can automatically detect the appropriate serializer based on the file ending, and saves in an atomic way by first saving to a temporary file and then moving it.

To define an encoding, an instance of :class:`.SerializerDispatch` is created, which has methods for saving and loading.

A default encoding which can handle common numpy and built-in types is also given.
"""

from ._version import __version__

from . import encoding
from ._iohelper import *
