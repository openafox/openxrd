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
from scipy.signal import argrelextrema

import matplotlib.pyplot as plt

import inspect


def find_peaks_2d(data, neigh_multi=0.01, neigh_add=5):

    neighborhood_size = int(data.size * neigh_multi) + neigh_add
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


def find_peaks_1d(data, neigh_multi=0.01, neigh_add=5):
    neighbors = int(data.size * neigh_multi) + neigh_add
    # print('neigh', neighbors)
    (x,) = argrelextrema(data, np.greater_equal, order=neighbors, mode='clip')
    return x.tolist()

def find_saddle_1d(data, neigh_multi=0.01, neigh_add=5):
    neighbors = int(data.size * neigh_multi) + neigh_add
    # print('neigh', neighbors)
    (x,) = argrelextrema(data, np.less_equal, order=neighbors)
    return x.tolist()


def fit_single(x, y, plot=False, model=PseudoVoigtModel, prefix=''):

    # run model
    mod = model(prefix=prefix)
    pars = mod.guess(y, x=x)
    out = mod.fit(y, pars, x=x)

    out = _out_addtion(x, y, out)

    #if plot:
    #    fig = plt.figure()
    #    out.plot(fig=fig)
    #    plt.show()

    if plot:
        args = {key: out.best_values[prefix + key] for key in
                inspect.getargspec(mod.func)[0] if key is not 'x'}
        plt.plot(x, model().func(x, **args), 'g--')
        plt.plot(x, y, 'b')
        plt.plot(x, out.init_fit, 'k--')
        plt.plot(x, out.best_fit, 'r-')
        # plt.savefig('../doc/_images/models_nistgauss2.png')
        plt.show()
    return out


def _out_addtion(x, y, out):
    """make new version of out from lmfit"""
    # Get fit values - could also use out.params[param].value
    # http://lmfit.github.io/lmfit-py/builtin_models.html#pseudovoigtmodel
    out.report = {}
    for key in out.params:
        #print(key)
        out.report[key] = out.params[key].value
        if out.params[key].stderr is not None:
            out.report[key+'_error'] = out.params[key].stderr
        else:
            out.report[key+'_error'] = 'N/A'

    nm = out.model._reprstring(long=True)
    nm = nm.replace("Model", "").replace("(", "").replace(")", "").split('+')
    out.report['mods'] = nm
    out.report['r^2'] = calc_r_sqd(y, out.best_fit)
    out.report.update({'height_obs':   np.max(y),
                       'y_min':        y.min(),
                       'mid_obs':      x[np.argmax(y)],
                        })
    return out


def calc_r_sqd(y, fit):
    """ Calculate the R^2 for a curve fit. """
    # https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html?s_tid=gn_loc_drop
    # https://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i
    # http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
    # http://kitchingroup.cheme.cmu.edu/blog/2013/02/18/Nonlinear-curve-fitting-with-confidence-intervals/
    SSe = np.sum((y - fit)**2.0)
    SSt = np.sum((y - np.mean(y))**2.0)
    return 1.0 - (SSe/SSt)


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
        rets = get_fit_all_1d(line, x_axis, row[1], plot)
        out[i].append(rets[0].params['fwhm'])
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

        out = get_fit(x_axis[left:right], line[left:right], plot)
        rets.append(out)

    return rets


def _set_bounds(x, y,  x_min, x_max, mids=None, num=1):
    """Set bounds if they are specified.
    Add mid bounds by finding:
    Saddles - mids='sad'
    Maximums - mids='max'
    None - mids=None
    """
    x_ = []
    # set x min if secified
    if x_min:
        x_.append(np.abs(x-x_min).argmin())
    else:
        x_.append(0)
    # set x max if secified
    if x_max:
        x_.append(np.abs(x-x_max).argmin())
    else:
        x_.append(len(x))

    # get mid points
    mid = []
    if mids in 'saddle':
        mid = find_saddle_1d(y)
    elif mids in 'maximum':
        mid = find_peaks_1d(y, 0.15)
    elif isinstance(mids, list):
        mid = mids
    # Remove edges from mid to prevent 0 array len errors
    mid = [i for i in mid if (abs(i-x_[0])>10 and abs(i-x_[-1])>10)]
    # print(mid)
    # Only and make sure to include needed number of mid points
    while num-1:
        val = sorted(mid)[(len(mid)-1)//2]
        x_.append(val)
        mid.remove(val)
        num += -1
        if not len(mid):
            mid = [x_[0]+(x_[0]+x_[-1])//2]
    x_.sort()
    # print(x_)
    return x_


def fit_multipeak(x, y,  name,
                  models=[PseudoVoigtModel, PseudoVoigtModel],
                  background_mod=None,
                  x_min=None, x_max=None, mids='max',
                  plot=False):

    x_ = _set_bounds(x, y, x_min, x_max, mids=mids, num=len(models))

    mod = []
    # Get the fits
    for i, model in enumerate(models):
        mod.append(model(prefix='mod_%d_' % i))
        #par = mod.guess(y[x1 + i * step: x1 + (i + 1) * step],
        #                x=x[x1 + i * step: x1 + (i + 1) * step])
        par = mod[i].guess(y[x_[i]:x_[i+1]], x=x[x_[i]:x_[i+1]])
        #par = mod.guess(y, x=x)
        par['mod_%d_amplitude' % i].set(min=1.e-13)
        if i < 1:
            mods = mod[i]
            pars = par
        else:
            mods += mod[i]
            pars += par

    if background_mod:
        i += 1
        if 'Poly' in background_mod.__name__:
            mod.append(background_mod(5, prefix='mod_%d_' % i))
        else:
            mod.append(background_mod(prefix='mod_%d_' % i))
        par = mod[i].guess(y, x=x)
        mods += mod[i]
        pars += par

    out = mods.fit(y[x_[0]:x_[-1]], pars, x=x[x_[0]:x_[-1]], nan_policy='omit')
    # save mods list
    out = _out_addtion(x, y, out)
    out.report['mod'] = mod
    if background_mod:
        if 'Poly' in background_mod.__name__:
            # account for other poly cases
            out.report['mod_%d_c6' % i] = 0
            out.report['mod_%d_c7' % i] = 0

    if plot:
        for i, model in enumerate(mod):
            args = {key: out.report['mod_%d_%s' % (i, key)] for key in
                    inspect.getargspec(model.func)[0] if key is not 'x'}
            plt.plot(x[x_[0]:x_[-1]], model.func(x[x_[0]:x_[-1]], **args),
                     'g--')
        plt.plot(x[x_[0]:x_[-1]], y[x_[0]:x_[-1]], 'b')
        plt.plot(x[x_[0]:x_[-1]], out.init_fit, 'k--')
        plt.plot(x[x_[0]:x_[-1]], out.best_fit, 'r-')
        # plt.savefig('../doc/_images/models_nistgauss2.png')
        plt.show()

    return out


def fits_to_csv(fits, keys, extra_data, name, savename):
    """Create CSV from fit report"""
    for i, out in enumerate(fits):
        k = 0
        out.report.update(extra_data[i])
        table = [keys]
        for j, modnm in enumerate(out.report['mods']):
            mod = 'mod_%d_' % j
            table.append([])
            k += 1
            for a, key in enumerate(keys):
                if key in 'name':
                    table[k].append(name)
                elif key in 'model':
                    table[k].append(modnm)
                elif mod+key in out.report:
                    table[k].append(out.report[mod+key])
                elif key in out.report:
                    table[k].append(out.report[key])
                elif '2d_error' in key:
                    col = chr(62+a)
                    col2 = chr(62+a+1)
                    table[k].append('=0.192575*SIN(%s2*PI()/360)*1'
                                    '/SIN(%s2*PI()/360)^3*%s2*PI()/180'
                                    % (col, col, col2))
                elif '2d' in key:
                    col = chr(63+a)
                    table[k].append('=1.540598/(SIN(%s2*PI()/360))' % col)
                else:
                    table[k].append('N/A')

        if os.path.exists(savename+'_newfits.csv'):
            table = table[1:]
        with open(savename+'_newfits.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerows(table)


def fits_to_csv_multitype(x, y,  name, savename, models=[PseudoVoigtModel],
                          x_min=None, x_max=None, psi=False,
                          extrahead=[], extra=[],
                          plot=False, plot_all=False, print_out=False):

    if len(extrahead) != len(extra):
        raise  Exception('extrahead and extra must be of same length')

    x1, x2 = _set_bounds(x, y, x_min, x_max)

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
        table[i] += [fits[j]['center'], fits[j]['center_error']]
        if not psi:
            table[i] += ['=1.540598/(SIN(%s2*PI()/360))' % col,
                         '=0.192575*SIN(%s2*PI()/360)*1/SIN(%s2*PI()/360)^3*%s2*PI()/180'
                         % (col, col, col2)]
    for j, mod_nm in enumerate(mod_nms):
        table[i] += [fits[j]['height'], fits[j]['fwhm'], fits[j]['amplitude']]
    for j, mod_nm in enumerate(mod_nms):
        table[i] += [fits[j]['sigma_error'], fits[j]['amplitude_error'],
                     fits[0]['r^2']]

    with open(savename+'fits.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerows(table)


def fit_data_to_csv(x, y,  name, savename,
                     x_min=None, x_max=None,
                     model=PseudoVoigtModel, plot=True):
    x1, x2 = _set_bounds(x, y, x_min, x_max)

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
    column = _check_array(column)
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
        if isinstance(column[0], list):
            table = column
        else:
            table = zip(column)
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(table)

def num_cols(array):
    if isinstance(array[0], list):
        return len(array[0])
    return 1

def _check_array(array):
    """Check if python list"""
    if not isinstance(array, list):
        if isinstance(array, np.ndarray):
            array = array.tolist()
        else:
            raise TypeError("input must be python list or numpy.ndarray")
    return array

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
