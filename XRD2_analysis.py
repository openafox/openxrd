#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
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
import inspect
from openxrd import get_datafiles
from openxrd import find_peaks_2d
from openxrd import find_peaks_1d
from openxrd import fit_single
from openxrd import get_fit_all_1d
from openxrd import get_fit_all_2d
from openxrd import fits_to_csv_multitype
from openxrd import fit_data_to_csv
from openxrd import BrukerData
from openxrd import fit_multipeak
from openxrd import fits_to_csv
from openxrd import csv_append_col

from datasetmetta import get_name_data

import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1 import make_axes_locatable

import csv
from lmfit import models
from lmfit import lineshapes

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm


def get_files():
    location = '/Users/towel/_The_Universe/_Materials_Engr/_Mat_Systems/_BNT_BKT/_CSD/_Data/EAPSI'
    # location = '.'
    files = get_datafiles(['*.raw'], location)
    return files


def merge_data(files):
    # Merge files into 1 for analysis
    data_list = []
    file_list = []
    for i, f in enumerate(sorted(files)):
        datafile = str(f)
        # Maybe add output that says joined bla
        # print(os.path.basename(datafile)[:-11])
        if len(data_list) < 1:
            data_list.append(BrukerData(datafile))
        else:
            data = BrukerData(datafile)
            diff = data.x[0] - data_list[-1].x[-1]
            step = data.rngs[0].metta['step_size']
            if (np.array_equal(data_list[-1].y, data.y) and
                    step*-2 < diff < 2*step):  # are ys same and xs overlap?
                data_list[-1] = data_list[-1] + data
            else:
                data_list.append(data)
        file_list.append(datafile)
    return data_list, file_list

def colorbar(mappable):
    #http://joseph-long.com/writing/colorbars/
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def plot_heatmap(data, title=None, mini=5, maxi=1e3, xy=None, plotpeaks=None):
    # colors
    # https://matplotlib.org/users/colormaps.html
    # https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.pcolor.html
    axlabels = {'family':               'sans',
                'color':                'black',
                'weight':               'normal',
                'size':                 20,
                }
    titles = {'family':                 'sans',
              'color':                  'black',
              'weight':                 'normal',
              'size':                   24,
              }
    labels = {'family':                 'sans',
              'fontname':               'DejaVu Sans',
              'color':                  '#66ff33',
              'weight':                 'normal',
              'size':                   14,
              'verticalalignment':      'center',
              'horizontalalignment':    'right'
              }

    fig = plt.figure(figsize=(10, 6), dpi=100)
    ax = fig.add_subplot(111)

    if title:
        ax.set_title(title, fontdict=titles)
    plot = ax.pcolormesh(data.x, data.y, data.smap, vmin=mini, vmax=maxi,
                         cmap='viridis')  # alpha=0.8)
    # plt.pcolor(x, y, data, norm=LogNorm(vmin=data.min()+5,
    #            vmax=data.max(), cmap='viridis') #alpha=0.8)
    ax.set_xlabel('2\u03b8[\u00b0]', fontdict=titles)
    ax.set_ylabel(u'\u03A8[\u00b0]', fontdict=titles)
    if xy is not None:
        points = ax.plot(xy[:, 1], xy[:, 0], 'ro', markersize=1)
    cbar = colorbar(plot)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["", ""])
    cbar.set_label("Intensity [arbitrary units]", fontdict=titles)
    # Set axis tick labels font
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        for prop in axlabels:
            getattr(label, 'set_' + prop)(axlabels[prop])
    if plotpeaks:
        # fastest?
        # https://softwarerecs.stackexchange.com/questions/7463/fastest-python-library-to-read-a-csv-file
        with open(plotpeaks, 'r') as f:
            peaks = csv.reader(row for row in f if not
                               row.startswith('#'))
            for peak in peaks:
                if (data.x.min() < float(peak[0]) < data.x.max() and
                        data.y.min() < float(peak[1]) < data.y.max()):
                    txt = ax.text(float(peak[0]), float(peak[1]), peak[2],
                                  fontdict=labels)
                    txt.set_path_effects([path_effects.Stroke(linewidth=1,
                                         foreground='black'),
                                         path_effects.Normal()])
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    """ The Following code is broken into sections that can be turned on or
    off depending on the desired analysis. Just change if to True or False.
    """

    # Get data
    files = get_files()
    data_list, file_list = merge_data(files)
    print(len(data_list))

    for i, data in enumerate(data_list):
        # Plot heat maps
        if False:
            xy_raw = find_peaks_2d(data.smap)
            # rescale xy peaks
            xy = np.asarray([data.get_real_xy(row[1], row[0])
                             for row in xy_raw])
            xy = np.roll(xy, 1, axis=1)  # quick fix. need to do properly

            mini = data.smap.min()
            ## make this base on fit of distribution??
            maxi = data.smap.max()*0.005
            peakfile = os.path.join(os.path.dirname(__file__),
                                    'BNKT_peaks.csv')
            plot_heatmap(data, maxi=maxi)
            # plot_heatmap(data, os.path.basename(file_list[i])[:19], maxi=maxi,
                         # xy=xy, plotpeaks=peakfile)
            plot_heatmap(data, maxi=maxi, plotpeaks=peakfile)

        # Do fits of all Automagically and save in CSV
        if False:
            # fit all 2th lines
            out_2th = get_fit_all_2d(data.smap, xy_raw, data.x, data.y,
                                     plot=False)
            # fit all psi lines
            smapT = data.smap.copy().T
            xy_raw = np.roll(xy_raw, 1, axis=1)
            out_psi = get_fit_all_2d(smapT, xy_raw, data.y, data.x, plot=False)
            # create table and write to csv
            table = []
            for i, row in enumerate(xy):
                table.append([row[0], row[1], out_2th[i][2], out_psi[i][2]])
            with open(datafile[:-11] + '.csv', 'wb') as f:
                writer = csv.writer(f)
                writer.writerows(table)

        # Manually fit a specific ranges
        if True:
            # insitue data
            if True:
                pos = []
                peaks = []
                # if specific to insitue flat measurement
                if (data.y[0] < 0 < data.y[-1] and
                        data.x[0] < 46 < data.x[-1]):
                    peaks.append('200')
                    pos.append(data.get_index_xy(46.79, -10))
                    pos.append(data.get_index_xy(41, 10))
                    pos.append(data.get_index_xy(53, 0))
                    """
                    peaks.append('110')
                    pos.append(data.get_index_xy(32.5, -10))
                    pos.append(data.get_index_xy(30, 10))
                    pos.append(data.get_index_xy(34, 0))
                    peaks.append('Pt111')
                    pos.append(data.get_index_xy(39.76, -10))
                    pos.append(data.get_index_xy(35, 10))
                    pos.append(data.get_index_xy(45, 0))
                    """
                # if specific to insitue tilted measurement
                elif (data.y[0] < 45 < data.y[-1] and
                      data.x[0] < 32 < data.x[-1]):
                    """
                    peaks.append('110_100')
                    pos.append(data.get_index_xy(32.56, 30))
                    pos.append(data.get_index_xy(30, 60))
                    pos.append(data.get_index_xy(35, 45))
                    peaks.append('111_100')
                    pos.append(data.get_index_xy(40.12, 50))
                    pos.append(data.get_index_xy(38, 60))
                    pos.append(data.get_index_xy(40.5, 54.74))
                    #pos.append(data.get_index_xy(35, 60))
                    #pos.append(data.get_index_xy(45, 54.74))
                    peaks.append('Pt200_111')
                    pos.append(data.get_index_xy(46.24, 45))
                    pos.append(data.get_index_xy(42, 65))
                    pos.append(data.get_index_xy(50, 54.74))
                    """
            # Si
            if False:
                pos = []
                peaks = []
                # if specific to insitue flat measurement
                peaks.append('Si_400')
                pos.append(data.get_index_xy(69.132, -10))
                pos.append(data.get_index_xy(68, 10))
                pos.append(data.get_index_xy(72, 0))
            # full map data
            if False:
                peakfile = os.path.join(os.path.dirname(__file__),
                                        'BNKT_peaks.csv')
                ## make peak file auto creat arrays as below
                peak = '100'
                pos.append(data.get_index_xy(22.93, -10))
                pos.append(data.get_index_xy(17, 10))
                pos.append(data.get_index_xy(29, 0))
                peak = '110'
                pos.append(data.get_index_xy(32.56, -10))
                pos.append(data.get_index_xy(30, 80))
                pos.append(data.get_index_xy(35, 45))

            smapT = data.smap.copy().T
            lines = []
            name = os.path.basename(file_list[i])[:23]
            basename = os.path.dirname(file_list[i])
            directory = os.path.abspath(os.path.join(basename, os.pardir))
            print(basename)
            print(name)
            for i in range(len(peaks)):
                print(peaks[i])
                # should we use an integrated area or a slice of the data
                integrate = True
                # which models to use
                # mods = [models.Pearson7Model, models.VoigtModel,
                #        models.PseudoVoigtModel]
                mods = [models.PseudoVoigtModel]
                # Psi
                if False:
                    x = data.y[pos[i*3][1]:pos[i*3+1][1]]
                    if integrate:
                        y = data.integrate_2d([pos[i*3+1][0], pos[i*3][1],
                                               pos[i*3+2][0], pos[i*3+1][1]],
                                               'y')
                    else:
                        y = smapT[pos[i*3][0], pos[i*3][1]:pos[i*3+1][1]]

                    # Do fit
                    savename = os.path.join(directory, '%s_psi' % peaks[i])
                    if True:
                        fits_to_csv_multitype(
                                    x, y, name, savename,  mods,
                                    psi=True,
                                    extrahead=['comp', 'thick', 'num', 'volt'],
                                    extra=get_name_data(name),
                                    plot=False, plot_all=False,
                                    print_out=False)

                    # Fit data to csv
                    if False:
                        fit_data_to_csv(x, y, name, savename, plot=False)

                csvheads = ['name', 'model',  'mid_obs', 'height_obs', 'comp',
                            'thick', 'num', 'volt', 'center','center_error',
                            '2d', '2d_error', 'height', 'fwhm', 'sigma',
                            'amplitude', 'sigma_error', 'amplitude_error',
                            'r^2', 'decay', 'slope', 'intercept']
                xdata = get_name_data(name)
                extradata = {'comp': xdata[0], 'thick': xdata[1],
                             'num': xdata[2], 'volt': xdata[3]}
                # 2th
                if True:
                    x = data.x[pos[i*3+1][0]:pos[i*3+2][0]]
                    if integrate:
                        y = data.integrate_2d([pos[i*3+1][0], pos[i*3][1],
                                               pos[i*3+2][0], pos[i*3+1][1]],
                                               'x')
                    else:
                        y = data.smap[pos[i*3+2][1], pos[i*3+1][0]:pos[i*3+2][0]]
                    # Do Fit
                    savename = os.path.join(directory, '%s_2th' % peaks[i])
                    if False:
                        fits_to_csv_multitype(
                                    x, y, name, savename,  mods,
                                    psi=False,
                                    extrahead=['comp', 'thick', 'num', 'volt'],
                                    extra=get_name_data(name),
                                    plot=True, plot_all=False,
                                    print_out=False)
                    # Fit data to csv
                    if False:
                        fit_data_to_csv(x, y, name, savename, plot=False)
                    # fit as multiple peaks
                    if True:
                    # compair these
                        out = []
                        rsqd = []
                        """
                        print('single')
                        out.append(fit_single(x, y, plot=False))
                        print('single poly')
                        out.append(fit_multipeak(x, y, name,
                                            models=[models.PseudoVoigtModel],
                                            background_mod=models.PolynomialModel,
                                            plot=False))
                        """
                        print('double poly')
                        out.append(fit_multipeak(x, y, name,
                                            models=[models.PseudoVoigtModel,
                                                    models.PseudoVoigtModel],
                                            background_mod=models.PolynomialModel,
                                            plot=False))
                        #fits_to_csv([o for o in out], csvheads,
                        #            [extradata, extradata, extradata],
                        #            name, savename)
                        # Add data to csv
                        fn = os.path.join(directory, 'fitdata')
                        csv_append_col(fn, [name + '_x'] + x.tolist())
                        csv_append_col(fn, [name + '_y'] + y.tolist())
                        csv_append_col(fn, [name + '_fit'] +
                                       out[0].best_fit.tolist())
                        for i, model in enumerate(out[0].report['mod']):
                            args = {key: out[0].report['mod_%d_%s' % (i, key)]
                                    for key in inspect.getargspec(model.func)[0]
                                    if key is not 'x'}
                            csv_append_col(fn, [name + '_f%d' %i] +
                                           model.func(x, **args).tolist())


    """
    ## Plot in difrent ways ######################################
    # surface #################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('2th')
    ax.set_zlabel('counts')
    ax.set_ylabel('psi')
    #ax.invert_zaxis()
    X, Y = np.meshgrid(data.x[x2:x3], data.y[y1:y2])
    Z = data.smap[y1:y2, x2:x3]
    ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    #ax.plot_wireframe(X, Y, Z, alpha=0.3, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='z', offset=0, cmap=cm.coolwarm)
    cset = ax.contourf(X, Y, Z, zdir='x', offset=39, cmap=cm.coolwarm)
    cset = ax.contourf(X, Y, Z, zdir='y', offset=5, cmap=cm.coolwarm)
    plt.show()


    # make 3d line plot of small range ##################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('2th')
    ax.set_zlabel('counts')
    ax.set_ylabel('psi')
    #ax.invert_zaxis()
    ys = data.y[y1:y2]
    for i, line in enumerate(data.smap[y1:y2, x2:x3]):
        y = [ys[i]] * len(line)
        ax.plot(data.x[x2:x3], y, line)
    #ax.view_init(120, 260)
    #plt.draw()
    plt.show()
    # The other way
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('2th')
    ax.set_zlabel('counts')
    ax.set_ylabel('psi')
    #ax.invert_zaxis()
    xs = data.x[x2:x3]
    for i, line in enumerate(smapT[x2:x3, y1:y2]):
        x = [xs[i]] * len(line)
        ax.plot(x, data.y[y1:y2], line)
    #ax.view_init(120, 260)
    #plt.draw()
    plt.show()
    """
    print('DONE')
