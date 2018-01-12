# -*- coding: utf-8 -*-
"""OpenXRD

This Package provides tools necessary to open many X-ray diffraction file
formats in an easilly manipulatible format.

Example:
    this will be an example::

        $ python example.py

Supported file types:
    Bruker RAW

Todo:
    * Solidify class structure
    * Finish Bruker to open all types including Gfrm
    * Add more file types
"""
# Copyright 2018 Austin Fox (openafox.com)
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

# Python 3 compatibility
from __future__ import (print_function, unicode_literals)
# #######################

from .databruker import BrukerData
from .analysistools import *

