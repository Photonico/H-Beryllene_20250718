#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    08.03.2016 10:21:06 CET
# File:    _save_load.py
"""Free functions for saving and loading objects with a given encoding."""

import os
import json
import pickle
import tempfile
from collections import namedtuple

import msgpack

from fsc.export import export

SerializerSpecs = namedtuple('SerializerSpecs',
                             ['binary', 'encode_kwargs', 'decode_kwargs'])


@export
class SerializerDispatch:
    """
    Defines functions for saving and loading for a given encoding.

    :param encoding: Object with ``encode`` and ``decode`` members, to be passed as the ``default`` and ``object_hook`` parameters, respectively (for JSON and msgpack, see e.g. :py:func:`json.dump` and :py:func:`json.load`).

    :param exclude: List of excluded serializers.
    :type exclude: list


    Basic usage:

    .. code :: python

        import fsc.iohelper

        IO_HANDLER = fsc.iohelper.SerializerDispatch(fsc.iohelper.encoding.default)
        IO_HANDLER.save([1, 2, 3], 'filename.json')
        x = IO_HANDLER.load('filename.json')
    """
    def __init__(self, encoding, exclude=[]):
        self.ext_mapping = {
            k.lower(): v
            for k, v in [('p', pickle), ('pickle',
                                         pickle), ('msgpack',
                                                   msgpack), ('json', json)]
            if v not in exclude
        }
        self.serializer_specs = dict([
            (json,
             SerializerSpecs(binary=False,
                             encode_kwargs=dict(default=encoding.encode),
                             decode_kwargs=dict(object_hook=encoding.decode))),
            (msgpack,
             SerializerSpecs(binary=True,
                             encode_kwargs=dict(default=encoding.encode),
                             decode_kwargs=dict(object_hook=encoding.decode))),
            (pickle,
             SerializerSpecs(binary=True,
                             encode_kwargs=dict(protocol=4),
                             decode_kwargs=dict()))
        ])
        for excluded_serializer in exclude:
            self.serializer_specs.pop(excluded_serializer, None)

    def _get_serializer(self, file_path, use_default=False):
        """Tries to determine the correct serializer from the file extension. If none can be determined, falls back to the default (:py:mod:`json`) if ``use_default`` is true, otherwise it throws an error."""
        _, file_ext = os.path.splitext(file_path)
        try:
            return self.ext_mapping[file_ext.lower().lstrip('.')]
        except KeyError:
            if use_default:
                return json
            else:
                raise ValueError(
                    "Could not guess serializer from file ending '{}'".format(
                        file_ext))

    def save(self, obj, file_path, serializer='auto'):
        """Saves an object to the file given in ``file_path``. The saving is made atomic (on systems where :py:func:`os.replace` is atomic) by first creating a temporary file and then moving to the ``file_path``.

        :param obj:         Object to be saved.

        :param file_path:   Path to the file.
        :type file_path:    str

        :param serializer:  The serializer to be used. Valid options are :py:mod:`msgpack` :py:mod:`json` and :py:mod:`pickle`. By default, the serializer is determined from the file extension. If this does not work, :py:mod:`json` is used to avoid data loss.
        :type serializer:   module
        """
        if serializer == 'auto':
            serializer = self._get_serializer(file_path, use_default=True)
        specs = self.serializer_specs[serializer]
        dirname = os.path.dirname(os.path.abspath(file_path))
        if not os.path.isdir(dirname):
            raise ValueError('Directory {} does not exist'.format(dirname))
        try:
            with tempfile.NamedTemporaryFile(
                    dir=dirname, delete=False,
                    mode='wb' if specs.binary else 'w') as f:
                tmp_path = f.name
                serializer.dump(obj, f, **specs.encode_kwargs)
            # closing necessary on Windows
            os.replace(tmp_path, file_path)
        except Exception as e:
            os.remove(tmp_path)
            raise e

    def load(self, file_path, serializer='auto'):
        """Loads the object that was saved to ``file_path``.

        :param file_path:   Path to the file.
        :type file_path:    str

        :param serializer:  The serializer which should be used to load the result. By default, is deduced from the file extension. If no serializer is given and it cannot be deduced from the file ending, a :py:class:`ValueError` is raised, to avoid loading corrupted data.
        :type serializer:   module
        """
        if serializer == 'auto':
            serializer = self._get_serializer(file_path)
        specs = self.serializer_specs[serializer]
        with open(file_path, 'rb' if specs.binary else 'r') as f:
            return serializer.load(f, **specs.decode_kwargs)
