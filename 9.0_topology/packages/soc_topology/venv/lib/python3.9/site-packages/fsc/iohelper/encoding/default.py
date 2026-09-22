#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    15.04.2016 10:21:56 CEST
# File:    default.py

"""
This is an example implementation of an encoding, to be used as input for :class:`.SerializerDispatch`.
"""

import numbers
import contextlib
from functools import singledispatch
from collections.abc import Iterable

import numpy as np
from fsc.export import export

@export
@singledispatch
def encode(obj):
    """
    Encodes numpy and built-in types into JSON / msgpack - compatible types.
    """
    raise TypeError('cannot JSONify {} object {}'.format(type(obj), obj))

@encode.register(np.bool_)
def _(obj):
    return bool(obj)

@encode.register(numbers.Real)
def _(obj):
    return float(obj)

@encode.register(numbers.Complex)
def _(obj):
    return dict(__complex__=True, real=encode(obj.real), imag=encode(obj.imag))

@encode.register(Iterable)
def _(obj):
    return list(obj)

#-----------------------------------------------------------------------#

@export
@singledispatch
def decode(obj):
    """
    Decodes JSON / msgpack objects into the corresponding built-in types.
    """
    return obj

def decode_complex(obj):
    return complex(obj['real'], obj['imag'])

@decode.register(dict)
def _(obj):
    with contextlib.suppress(AttributeError):
        obj = {k.decode('utf-8'): v for k, v in obj.items()}
    special_markers = [key for key in obj.keys() if key.startswith('__')]
    if len(special_markers) == 1:
        name = special_markers[0].strip('__')
        return globals()['decode_' + name](obj)
    else:
        return obj
