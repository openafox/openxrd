#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Opens Bruker RAW V3, 2D and 3D Data and stores in nice format - based on
"General Area Detector Diffraction System (GADDS) Version 4.1.xx"
Apendix B. and https://github.com/wojdyr/xylib"""
# Copyright 2015 Austin Fox
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

import sys, os
import numpy as np
import struct


class BrukerHeader(object):
    """Bruker Raw Header Dictonary"""

    def __init__(self):
        # For each metta data set keywords are contain (1) the actual
        # value, (2) a user-suitable label for the item, (3) the data type
        # and (4) the data position in the group:
        # ordering index:

        self.attrs = {
            'version':       [None, 'Version',                     4,   0],
            'head_1':        [None, 'head 1',                      4,   4],
            'file_status':   [None, 'File Status',              '<I',   8],
            'range_cnt':     [None, 'Range Count',              '<I',  12],
            'm_date':        [None, 'Measure Date',               10,  16],
            'm_time':        [None, 'Measure Time',               10,  26],
            'user':          [None, 'User',                       72,  36],
            'site':          [None, 'Site',                      218, 108],
            'sample_id':     [None, 'Sample ID',                  60, 326],
            'comment':       [None, 'Comment',                   160, 386],
            'head_2':        [None, 'head 2',                      2, 546],
            'c_goni':        [None, 'Goniometer Model',         '<I', 548],
            'c_goni_s':      [None, 'Goniometer Stage',         '<I', 552],
            'c_samp_l':      [None, 'Sample Changer',           '<I', 556],
            'c_goni_c':      [None, 'Goniometer Controler',     '<I', 560],
            'c_goni_r':      [None, '(R4) goniometer radius',   '<f', 564],
            'fix_divr':      [None, '(R4) fixed divergence',    '<f', 568],
            'fix_samp':      [None, '(R4) fixed sample slit',   '<f', 572],
            'prim_ss':       [None, 'primary Soller slit',      '<I', 576],
            'prim_mon':      [None, 'primary monochromator',    '<I', 580],
            'fix_anti':      [None, '(R4) fixed antiscatter',   '<f', 584],
            'fix_detc':      [None, '(R4) fixed detector slit', '<f', 588],
            'sec_ss':        [None, 'secondary Soller slit',    '<f', 592],
            'fix_tf':        [None, 'fixed thin film attach',   '<I', 596],
            'beta_f':        [None, 'beta filter',                 4, 600],
            'sec_mon':       [None, 'secondary monochromator',  '<f', 604],
            'anode':         [None, 'Anode Material',              4, 608],
            'head_3':        [None, 'head 3',                      4, 612],
            'alpha_ave':     [None, 'Alpha Average',            '<d', 616],
            'alpha_1':       [None, 'Alpha 1',                  '<d', 624],
            'alpha_2':       [None, 'Alpha 2',                  '<d', 632],
            'beta':          [None, 'Beta',                     '<d', 640],
            'alpha_ratio':   [None, 'Alpha_ratio',              '<d', 648],
            'unit_nm':       [None, '(C4) Unit Name',              4, 656],
            'int_beta_a1':   [None, 'Intensity Beta:a1',           4, 660],
            'mea_time':      [None, 'Measurement Time',         '<f', 664],
            'head_4':        [None, 'head 4',                     43, 668],
            'hard_dep':      [None, 'hard_dep',                    1, 711],
            }

    def __getitem__(self, key):
        if key in self.attrs:
            return self.attrs[key][0]
        return None

    def __len__(self):
        return len(self.attrs)

    def __setitem__(self, key, item):
        if key in self.attrs:
            self.attrs[key][0] = item
        else:
            self.attrs[key] = [item, key, len(self.attrs)]

    def __delitem__(self, key):
        if key in self.attrs:
            del self.attrs[key][0]

    def __iter__(self):
        self.index = -1
        return self

    def __next__(self):
        # Get items sorted in specified order:
        keys = sorted(self.attrs, key=lambda key: self.attrs[key][3])
        if self.index == len(keys) - 1:
            raise StopIteration
        self.index = self.index + 1
        key = keys[self.index]
        return key


class BrukerRangeHeader(BrukerHeader):
    """Bruker Raw Header Dictonary"""

    def __init__(self):
        # For each metta data set keywords are contain (1) the actual
        # value, (2) a user-suitable label for the item, (3) the data type
        # and (4) the data position in the group:
        # ordering index:

        self.attrs = {
            'header_len':    [None, 'Header Len',               '<I',   0],
            'steps':         [None, 'Steps',                    '<I',   4],
            'start_theta':   [None, 'Start Theta',              '<d',   8],
            'start_2th':     [0,    'Start 2Theta',             '<d',  16],
            'drive_chi':     [0,    'Chi Start',                '<d',  24],
            'drive_phi':     [None, 'Phi Start',                '<d',  32],
            'drive_x':       [None, 'X Start',                  '<d',  40],
            'drive_y':       [None, 'Y Start',                  '<d',  48],
            'drive_z':       [None, 'Z Start',                  '<d',  56],
            'ig_1':          [None, 'ig 1',                     '<Q',  64],
            'ig_2':          [None, 'ig 2',                        6,  72],
            'ig_2_1':        [None, 'ig 2_1',                   '<h',  78],
            'R8':            [None, '(R8) variable anitscat',   '<d',  80],
            'ig_3':          [None, 'ig 3',                        6,  88],
            'ig_3_1':        [None, 'ig 3_1',                   '<h',  94],
            'dec_code':      [None, 'Detector',                 '<I',  96],
            'hv':            [None, 'High Voltage',             '<f', 100],
            'amp_gain':      [None, 'Ampliphier Gain',          '<f', 104],
            'dis1_LL':       [None, 'Discriminator1 Lower Lev', '<f', 112],
            'ig_4':          [None, 'ig 4',                     '<I', 116],
            'ig_5':          [None, 'ig 5',                     '<d', 120],
            'ig_6':          [None, 'ig 6',                     '<f', 128],
            'ig_a':          [None, 'ig a',                     '<f', 132],
            'ig_b':          [None, 'ig b',                        5, 136],
            'ig_b_1':        [None, 'ig b_1',                      3, 141],
            'ig_c':          [None, 'Aux Axis 1 start',         '<d', 144],
            'ig_d':          [None, 'Aux Axis 2 start',         '<d', 152],
            'ig_e':          [None, 'Aux Axis 3 start',         '<d', 160],
            'ig_f':          [None, 'Scan Mode',                   4, 168],
            'ig_g':          [None, 'ig g',                     '<I', 172],
            'ig_h':          [None, 'ig h',                     '<I', 172],
            'step_size':     [None, 'Step Size',                '<d', 176],
            'ig_i':          [None, 'ig i',                     '<d', 184],
            'step_time':     [None, 'Time Per Step',            '<f', 192],
            'ig_j':          [None, 'Scan Type',                '<I', 196],
            'ig_k':          [None, 'Delay Time',               '<f', 200],
            'ig_l':          [None, 'ig l',                     '<I', 204],
            'rot_speed':     [None, 'Rotation Speed',           '<f', 208],
            'ig_m':          [None, 'ig m',                     '<f', 212],
            'ig_n':          [None, 'ig n',                     '<I', 216],
            'ig_o':          [None, 'ig o',                     '<I', 220],
            'gen_v':         [None, 'Generator Voltage',        '<I', 224],
            'gen_a':         [None, 'Generator Current',        '<I', 228],
            'ig_p':          [None, 'ig p',                     '<I', 232],
            'ig_q':          [None, 'ig q',                     '<I', 236],
            'lambda':        [None, 'Lambda',                   '<d', 240],
            'ig_r':          [None, 'ig r',                     '<I', 248],
            'ig_s':          [None, 'Len of each data in bits', '<I', 252],
            'sup_len':       [None, 'supplementary header len', '<I', 256],
            'ig_t':          [None, 'ig t',                     '<I', 260],
            'ig_u':          [None, 'ig u',                     '<I', 264],
            'ig_v':          [None, 'ig v',                     '<I', 268],
            'ig_w':          [None, 'ig w',                     '<I', 272],
            'ig_x':          [None, 'Reserved for expansion',   '<I', 280],
            }


class BrukerSupplementalHeader(BrukerHeader):
    """Bruker Raw Supplemental Header Dictonary"""

    def __init__(self):
        # For each metta data set keywords are contain (1) the actual
        # value, (2) a user-suitable label for the item, (3) the data type
        # and (4) the data position in the group:
        # ordering index:
        self.attrs = {
            'type':          [None, 'Record type',              '<I',   0],
            'length':        [None, 'record length',            '<I',   4],
            'reserved':      [None, 'reserved',                 '<I',   8],
            'int_start':     [None, 'integration range start',  '<f',  16],
            'int_end':       [None, 'integration range end',    '<f',  20],
            'chi_start':     [None, 'int range chi start',      '<f',  24],
            'chi_end':       [None, 'int range chi end',        '<f',  28],
            'norm':          [None, 'Normalization method',     '<I',  32],
            'prog':          [None, 'program name',               20,  36],
            'act_2th':       [None, 'act 2th',                  '<f',  56],
            'act_omega':     [None, 'act omega',                '<f',  60],
            'act_phi':       [None, 'act phi',                  '<f',  64],
            'act_psi':       [None, 'act psi',                  '<f',  68],
            }


class BrukerRange(object):

    def __init__(self):
        self.metta = {}
        self.supmetta = {}
        self.counts_data = []


class BrukerData(object):
    """Retrieves and stores Bruker XRD Data"""

    def __init__(self, filename=None):
        self.filename = filename
        self.rngs = []
        if filename:
            self.filecontent = self.get_data_from_file(self.filename)
            self.header = self.get_metta(BrukerHeader(), 0)
            pos = 712
            for i in range(self.header['range_cnt']):
                rng, pos = self.get_range(pos)
                self.add_range(rng)
            if self.rngs[0].supmetta['type'] == 200:  # Area map
                self.x = []
                self.y = []
                self.get_smap()
            else:
                raise Exception("not file from area detector, this is "
                                "currently not supported. Sorry")
        else:
            self.header = None
            self.x = []
            self.y = []
            self.smap = np.array([[]])

    def add_range(self, rng):
        self.rngs.append(rng)

    def get_data_from_file(self, filename):
        try:
            with open(filename, mode='rb') as f:  # b is important -> binary
                filecontent = f.read()
        except IOError as e:
            raise Exception("I/O error({0}): {1}".format(e.errno, e.strerror))

        if b"RAW1.01" not in filecontent[0:8]:
            raise Exception("invalid file type must be RAW1.01")
            #add other versions and gfrm
        return filecontent

    def get_range(self, pos):
        rng = BrukerRange()
        rng.metta = self.get_metta(BrukerRangeHeader(), pos)
        pos += rng.metta['header_len']
        rng.counts_data = []
        (typ, ) = struct.unpack('<I', self.filecontent[pos: pos+4])
        # Check type of supplemental
        if typ == 200:  # Area Detector Parameters
            rng.supmetta = self.get_metta(BrukerSupplementalHeader(), pos)
        elif typ == 190:  # offset assigned by EVA
            pass
        elif typ == 150:  # removed data for search
            pass
        elif typ == 140:  # comment
            pass
        elif typ == 130:  # QCI parameters (obsolete)
            pass
        elif typ == 120:  # OQM parameters
            pass
        elif typ == 110:  # PSD parameters
            pass
        elif typ == 100:  # oscillation parameters
            pass
        pos += rng.metta['sup_len']
        data_len = rng.metta['steps']
        for i in range(data_len):
            (ret,) = struct.unpack('<f', self.filecontent[pos+i*4: pos+i*4+4])
            rng.counts_data.append(ret)
        pos += data_len * 4
        return rng, pos

    def get_smap(self):
        y = []
        smap = []

        for i, rng in enumerate(self.rngs):
            phi = rng.supmetta['act_psi']
            chi_s = rng.supmetta['chi_start']
            chi_e = rng.supmetta['chi_end']
            y.append(90 - phi +
                     (90+chi_s-(chi_s-chi_e)/2))

            # print(phi, chi_s, chi_e, y[i])
            smap.append(rng.counts_data)
        # Sadly must handel silyness of bruker not keeping ranges same length
        # https://stackoverflow.com/questions/27890052/convert-and-pad-a-list-to-numpy-array
        lens = np.array([len(item) for item in smap])
        mask = lens[:, None] > np.arange(lens.max())
        # print(mask.shape)
        out = np.full(mask.shape, 0)
        out[mask] = np.concatenate(smap)

        userng = np.argmax(lens)
        twoth_0 = self.rngs[userng].metta['start_2th']
        twoth_s = self.rngs[userng].metta['step_size']
        twoth_e = twoth_0 + twoth_s*self.rngs[userng].metta['steps']
        twoth_len = lens.max()
        self.x = np.linspace(twoth_0, twoth_e, twoth_len)
        self.y = np.linspace(min(y), max(y), len(self.rngs))
        # convert smap to np.array
        self.smap = np.asarray(out)
        # print(self.x.shape, self.y.shape, self.smap.shape)

    def get_real_xy(self, x, y):
        y_out = self.y[int(y)]
        x_out = self.x[int(x)]
        return x_out, y_out

    def get_index_xy(self, x, y):
        """Assumes x and y are ordered arrays with len > 0"""
        y_out = np.abs(self.y-y).argmin()
        if len(self.y) == 1 and y > self.y[0]:
            y_out = 1
        x_out = np.abs(self.x-x).argmin()
        if len(self.x) == 1 and x > self.x[0]:
            x_out = 1
        return x_out, y_out

    def integrate_2d(self, area='all', axis='x'):
        """ area - (x1, y1, x2, y2)
            axis - axis to preserve"""
        if area == 'all':
            x1 = 0
            y1 = 0
            (y2, x2) = self.smap.shape
        elif len(area) != 4:
            raise Exception('area must be "all" or (x1, y1, x2, y2).')
        else:
            x1, y1, x2, y2 = area

        if axis == 'x':
            line = np.sum(self.smap[y1:y2, x1:x2], axis=0)  # 2th
        elif axis == 'y':
            line = np.sum(self.smap[y1:y2, x1:x2], axis=1)  # psi
        else:
            raise Exception('axis must be specified either "x" or "y".')
        return line

    def get_metta(self, mettaclass, start_pos):

        for key in mettaclass:
            pos = mettaclass.attrs[key][3]+start_pos
            typ = mettaclass.attrs[key][2]
            bits = 0
            if typ == '<h' or typ == '<H':
                bits = 2
            elif typ == '<f' or typ == '<I':
                bits = 4
            elif typ == '<d' or typ == '<Q' or typ == '<q':
                bits = 8
            elif isinstance(typ, int):
                (mettaclass[key],) = struct.unpack('%ds' % typ,
                            self.filecontent[pos: pos+typ])
                continue
            (mettaclass[key],) = struct.unpack(typ,
                            self.filecontent[pos: pos+bits])
        return mettaclass

    def __add__(self, other):
        try:
            if not np.array_equal(self.y, other.y):
                raise Exception('Must have same y scale')
        except:
            raise Exception('Type mismatch, must be BrukerData')
        ret = BrukerData()
        ret.y = self.y
        ret.smap = np.concatenate((self.smap, other.smap), axis=1)
        ret.x = np.concatenate((self.x, other.x), axis=0)
        ret.rngs = self.rngs + other.rngs
        ret.filename = self.filename + " & " + other.filename

        return ret

    # def __radd__(self, other):


if __name__ == '__main__':
    pass
