#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2018 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

# Python 3 compatibility
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (
         bytes, dict, int, list, object, range, str,
         ascii, chr, hex, input, next, oct, open,
         pow, round, super,
         filter, map, zip)
# #######################

def get_name_data(name):

    comp_thick = {"19_24": (2.5, 200),
                  "20_19": (2.5, 400),
                  "20_21": (2.5, 400),
                  "21_23": (2.5, 600),
                  "21_09": (5, 200),
                  "21_10": (5, 400),
                  "20_08": (5, 400),
                  "21_12": (5, 400),
                  "21_29": (5, 600),
                  "19_15": (10, 200),
                  "20_15": (10, 400),
                  "20_13": (10, 400),
                  "21_36": (10, 600),
                  "19_19": (0, 200),
                  "20_32": (0, 400),
                  "20_33": (0, 400),
                  "19_32": (0, 400),
                  "19_33": (0, 400),
                  "20_31": (0, 400),
                  "20_03": (0, 600)
                  }

    num = None
    volt = None
    if len(name) < 9:
        comp, thick = comp_thick[name[0:2] + '_' + name[2:4]]
        num = name[5:9]
    if len(name) > 8:
        comp, thick = comp_thick[name[8:13]]
    if len(name) > 16:
        num = name[17:19]
        volt = name[20:23]
    return [comp, thick, num, volt]


if __name__ == '__main__':
    pass
