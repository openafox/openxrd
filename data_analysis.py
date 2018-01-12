#!/usr/bin/env python
"""Data analysis functions for XRD2 """
# Copyright 2017 Austin Fox
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
import csv
from lmfit.models import PseudoVoigtModel
import lmfit
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from scipy.signal import argrelmax

import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog


def get_datafiles(supported_datafiles, location):
    """Qt file dialogue widget
    """
    types = ' '.join([row[0] for row in supported_datafiles])
    filetypes = 'Supported (' + types + ')'
    app = QApplication(sys.argv)
    widget = QWidget()
    files, _ = QFileDialog.getOpenFileNames(
                        widget,
                        'Program to run',
                        location,
                        filetypes + ';;All files (*.*)',
                        None,
                        QFileDialog.DontUseNativeDialog)
    return files


def find_peaks_2d(data):

    neighborhood_size = data.size * 0.0001
    threshold = np.average(data) * 2  # 1500
    # print(neighborhood_size, threshold)

    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    xy = np.array(ndimage.center_of_mass(data, labeled,
                                         range(1, num_objects+1)))
    return xy


def find_peaks_1d(data, multi=0.01, add=5):
    neighbors = int(data.size * multi) + add
    # print('neigh', neighbors)
    (x,) = argrelmax(data, order=neighbors)
    return x


def get_fit(x, y, plot=False, model=PseudoVoigtModel):

    ret = {}
    # first zero background
    ret['y_min'] = y.min()
    y = y-ret['y_min']
    # run model
    pars = model().guess(y, x=x)
    out = model().fit(y, pars, x=x)
    # Get fit values - could also use out.params[param].value
    # http://lmfit.github.io/lmfit-py/builtin_models.html#pseudovoigtmodel
    for key in out.params:
        ret[key] = out.params[key].value
    ret.update({'height_obs': np.max(y),
                'mid_obs': x[np.argmax(y)],
                'fit': out.best_fit,
                'full': out
                })
    pcov = out.covar
    perror = []
    i = 0
    for key in out.params.keys():
        if 'fwhm' not in key and 'height' not in key and 'gamma' not in key:
            if pcov is not None:
                ret[key[0:3]+'_error'] = np.absolute(pcov[i][i])**0.5
            else:
                ret[key[0:3]+'_error'] = 'N/A'
            i += 1

    # https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html?s_tid=gn_loc_drop
    SSe = np.sum((y - ret['fit'])**2.0)
    SSt = np.sum((y - np.mean(y))**2.0)
    ret['r^2'] = 1.0 - (SSe/SSt)
    # https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i
    # http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
    # http://kitchingroup.cheme.cmu.edu/blog/2013/02/18/Nonlinear-curve-fitting-with-confidence-intervals/
    if plot:
        fig = plt.figure()
        out.plot(fig=fig)
        plt.show()
    return ret


def div0(a, b):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


def get_fit_all_2d(datas, xy_raw, x_axis, y_axis, plot=False):
    j = 0
    group = []
    line = []
    out = xy_raw.tolist()
    for i, row in enumerate(xy_raw):
        oldline = line
        line = datas[int(row[0]), :]
        ret = get_fit_all_1d(line, x_axis, row[1], plot)
        out[i].append(ret[0]['fwhm'])
        """
            if xy_raw[j][0] == row[0]:  # get all on same line
                group.append([x_axis[left:right], ret['fit']])
                j = i
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.semilogy(x_axis, oldline, label=y_axis[int(row[0])])  # plot
                for test in group:
                    ax.semilogy(test[0], test[1], label='fit')
                j = i
                plt.show()
                group = []
                group.append([x_axis[left:right], ret['fit']])
        """
    return out


def get_fit_all_1d(line, x_axis, position=None, maxs=None, plot=False):
    if maxs is None:
        maxs = find_peaks_1d(line)
    print('maxs', maxs)
    rets = []
    if position is None:
        positions = maxs
    else:
        positions = [position]
    for pos in positions:
        m_i = (np.abs(maxs-int(pos))).argmin()
        if m_i > 0:
            # find min between desired pos and max before
            left = (int(maxs[m_i-1]) +
                    np.argmin(line[int(maxs[m_i-1]):int(pos)]))
        elif int(pos == 0):
            left = 0
        else:
            left = np.argmin(line[0:int(pos)])
        if m_i < len(maxs) - 1:
            right = (int(pos) +
                     np.argmin(line[int(pos):int(maxs[m_i+1])]))
        else:
            right = len(line) - 1

        if right - left < 4:
            right += 4 - (right - left)

        if plot:
            ret = get_fit(x_axis[left:right], line[left:right], True)
        else:
            ret = get_fit(x_axis[left:right], line[left:right], False)
        rets.append(ret)

    return rets


def fits_to_csv_multitype(x, y,  name, savename, models=[PseudoVoigtModel],
                          x_min=None, x_max=None, psi=False,
                          extrahead=[], extra=[],
                          plot=False, plot_all=False, print_out=False):

    if len(extrahead) != len(extra):
        raise  Exception('extrahead and extra must be of same length')

    # set limits if secified
    if x_min:
        x1 = np.abs(x-x_min).argmin()
    else:
        x1 = 0
    if x_max:
        x2 = np.abs(x-x_max).argmin()
    else:
        x2 = len(x) - 1

    fits = []
    mod_nms = []
    # Get the fits
    for model in models:
        mod_nms.append(model.__name__[0:-5])
        fits.append(get_fit(x[x1:x2], y[x1:x2], plot=plot, model=model))

    # Print all fits as a Table
    if print_out:
        print(' '*20, ''.join('{:^20s}'.format(s) for s in mod_nms))
        all_keys = [k for d in fits for k in d.keys()]
        keys = {k for k in all_keys if all_keys.count(k)==len(fits)}
        for key in keys:
            if key != 'full' and key != 'fit':
                print('{:^20s}'.format(key),
                      ''.join('{:^20f}'.format(fit[key]) for fit in fits))
    # plot all together
    if plot_all:
        fig, ax = plt.subplots()
        ax.plot(x[x1:x2], y[x1:x2]-min(y))
        for j, mod_nm in enumerate(mod_nms):
            ax.plot(x[x1:x2], fits[j]['fit'], label=mod_nm)
        ax.legend()
        # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", ncol=3)
        ax.set_ylabel('Intensity')
        if psi:
            ax.set_xlabel(u'\u03a8[\u00b0]')
        else:
            ax.set_xlabel(u'2\u03b8[\u00b0]')
        plt.show()

    # make table for csv (really a row unless adding header)
    table = []
    i = 0
    if not os.path.exists(savename+'fits.csv'):
        # Add headers
        table = [['name', 'mid_obs', 'height_obs'] + extrahead]
        for mod_nm in mod_nms:
            table[0] += [mod_nm + '_mid', mod_nm + 'mid_err']
            if not psi:
                table[0] += [mod_nm + '_2d',  mod_nm + '_2d_err']
        for mod_nm in mod_nms:
            table[0] += [mod_nm + '_height', mod_nm + '_fwhm', mod_nm + '_Area']
        for mod_nm in mod_nms:
            table[0] += [mod_nm + '_sig_err', mod_nm + '_A_err', mod_nm + '_R^2']

        i = 1

    # Add Data
    table.append([name, fits[0]['mid_obs'], fits[0]['height_obs']] + extra)
    for j, mod_nm in enumerate(mod_nms):
        col = chr(65+3+j*4+len(extra))
        col2 = chr(65+4+j*4+len(extra))
        table[i] += [fits[j]['center'], fits[j]['cen_error']]
        if not psi:
            table[i] += ['=1.540598/(SIN(%s2*PI()/360))' % col,
                         '=0.192575*SIN(%s2*PI()/360)*1/SIN(%s2*PI()/360)^3*%s2*PI()/180'
                         % (col, col, col2)]
    for j, mod_nm in enumerate(mod_nms):
        table[i] += [fits[j]['height'], fits[j]['fwhm'], fits[j]['amplitude']]
    for j, mod_nm in enumerate(mod_nms):
        table[i] += [fits[j]['sig_error'], fits[j]['amp_error'],
                     fits[0]['r^2']]

    with open(savename+'fits.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerows(table)


def fit_data_to_csv(x, y,  name, savename,
                     x_min=None, x_max=None,
                     model=PseudoVoigtModel, plot=True):

    if x_min:
        x1 = np.abs(x-x_min).argmin()
    else:
        x1 = 0
    if x_max:
        x2 = np.abs(x-x_max).argmin()
    else:
        x2 = len(x) - 1

    fit = get_fit(x[x1:x2], y[x1:x2], plot=plot, model=model)

    # Append column to csv
    if not os.path.exists(savename+'fitdata.csv'):
        table = [[x[x1:x2][j]] + [y[x1:x2][j]] + [row] for j, row in
                 enumerate(fit['fit'].tolist())]
        table.insert(0, ['angle', name, name + '_fit'])
        with open(savename+'fitdata.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(table)
    else:
        csv_append_col(savename+'fitdata', [name] + y[x1:x2].tolist())
        csv_append_col(savename+'fitdata', [name + '_fit'] +
                       fit['fit'].tolist())


def csv_append_col(filename, column, len_policy="ls"):
    """Append columns to a csv file.
    len_policy allows column(s) to add to be longer 'l', shorter 's', or
    require the exact lengths 'not l or s' as the table in the file.
    'ls' or 'sl' is also valid.
    """
    if not isinstance(column, list):
        if isinstance(column, np.ndarray):
            column = column.tolist()
        else:
            raise TypeError("column must be python list or numpy.ndarray")

    if not '.csv' in filename[-4:]:
        filename += '.csv'
    if os.path.exists(filename):
        with open(filename, 'r') as f_in:
            filetable = list(csv.reader(f_in))
            f_len = len(filetable)
            f_cols = num_cols(filetable)
            c_len = len(column)
            c_cols = num_cols(column)
            if f_len > c_len:
                if "s" in len_policy:
                    addit = [[""] * c_cols] * (f_len - c_len)
                    column = column + addit
                else:
                    raise ValueError("column (len=%d) can not be shorter then "
                            "the table in the file (len=%d), unless "
                            "len_policy='s'" % (c_len, f_len))

            if f_len < c_len:
                if "l" in len_policy:
                    addit = [[""] * f_cols] * (c_len - f_len)
                    cols = num_cols(filetable)
                    filetable += (addit)
                else:
                    raise ValueError("column (len=%d) can not be longer then "
                            "the table in the file (len=%d), unless "
                            "len_policy='l'" % (c_len, f_len))

            if c_cols == 1 and not isinstance(column[0], list):
                table = [row + [column[j]] for j, row in enumerate(filetable)]
            else:
                table = [row + column[j] for j, row in enumerate(filetable)]

    else:
        table = column
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(table)

def num_cols(array):
    if isinstance(array[0], list):
        return len(array[0])
    return 1


def fits_to_csv_autopeaks(x, y, savename):
    maxs = find_peaks_1d(y, 0.20)
    rets = get_fit_all_1d(y, x, maxs=maxs, plot=False)
    print('rets', len(rets))
    table = [['cent', 'fwhm', 'height_obs', 'height']]
    for fit in rets:
        table.append([fit['cent'], fit['fwhm'], fit['height_obs'],
                     fit['height']])

    with open(savename+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(table)


# Future adds deconvolving peaks
# http://kitchingroup.cheme.cmu.edu/blog/2013/01/29/Curve-fitting-to-get-overlapping-peak-areas/
# https://stackoverflow.com/questions/10143905/python-two-curve-gaussian-fitting-with-non-linear-least-squares
# and for ideas
# https://www.wavemetrics.com/products/igorpro/dataanalysis/peakanalysis/multipeakfitting.htm

# explination of how to calc errors - propegation
# https://chem.libretexts.org/Core/Analytical_Chemistry/Quantifying_Nature/Significant_Digits/Propagation_of_Error


if __name__ == '__main__':
    # datafile = get_datafiles([".txt"])

    # Test csv_append_col #
    exfiles = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                                           'examplefiles'))
    crnt_file = os.path.join(exfiles, "CRNT.csv")
    with open(crnt_file, 'r') as f_in:
            reader = csv.reader(f_in)
            col_len = len(list(reader))
    add_col = []
    add_col.append(list(range(0, col_len - 10)))
    add_col.append([[x] for x in range(0, col_len)])
    add_col.append(np.linspace(0, 10, col_len - 10).tolist())
    # add_col.append(np.linspace(0, 10, col_len - 10).tolist())
    for col in add_col:
        csv_append_col(crnt_file, col, '')
        print('try')
    #add_col = list(map(list, zip(*add_col)))
    #print(len(add_col), len(add_col[0]))
    #csv_append_col(crnt_file, add_col)
    #print('done')
    #########################
