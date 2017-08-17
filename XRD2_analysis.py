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
from data_analysis import get_datafiles
from data_analysis import find_peaks_2d
from data_analysis import find_peaks_1d
from data_analysis import get_fit
from data_analysis import get_fit_all_2d
from data_analysis import get_fit_all_1d
from data_analysis import fits_to_csv2
from bruker_data import BrukerData
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm


def get_files():
    location = '~/_The_Universe/_Materials_Engr/_Mat_Systems/_BNT_BKT/_CSD/'
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
                    step*-2 < diff < 2*step):
                data_list[-1] = data_list[-1] + data
            else:
                data_list.append(data)
        file_list.append(datafile)
    return data_list, file_list


files = get_files()
data_list, file_list = merge_data(files)
"""
# Find peaks
xy_raw = find_peaks_2d(data.smap)
# rescale xy peaks
xy = np.asarray([data.get_real_xy(row[1], row[0]) for row in xy_raw])
xy = np.roll(xy, 1, axis=1)  # quick fix. need to do properly
"""

peakfile = os.path.join(os.path.dirname(__file__), 'BNKT_peaks.csv')
print(len(data_list))

for i, data in enumerate(data_list):
    #""" Plot heat maps
    #data.plot_heatmap(name)
    mini=data.smap.min()
    maxi=data.smap.max()*.1
    data.plot_heatmap(os.path.basename(file_list[i])[:19], maxi=maxi,
                        plotpeaks=peakfile)
    #""

    """ Do fits of all Automagically and save in CSV
    # fit all 2th lines
    out_2th = get_fit_all_2d(data.smap, xy_raw, data.x, data.y, plot=False)
    # fit all psi lines
    smapT = data.smap.copy().T
    xy_raw = np.roll(xy_raw, 1, axis=1)
    out_psi = get_fit_all_2d(smapT, xy_raw, data.y, data.x, plot=False)
    table = []
    for i, row in enumerate(xy):
        table.append([row[0], row[1], out_2th[i][2], out_psi[i][2]])
    with open(datafile[:-11] +'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(table)
    """

    # """ Get and fit a range from multiple datafiles and save to csv
    # Manually fit a range
    """
    peak = '100'
    x1, y1 = data.get_index_xy(22.93, -10)
    x2, y2 = data.get_index_xy(17, 10)
    x3, y3 = data.get_index_xy(29, 0)
    """
    xs = [0]*9
    ys = [0]*9
    peaks = []
    if data.y[0] < 0 < data.y[-1] and data.x[0] < 46 < data.x[-1]:
        peaks.append('200')
        xs[0], ys[0] = data.get_index_xy(46.79, -10)
        xs[1], ys[1] = data.get_index_xy(41, 10)
        xs[2], ys[2] = data.get_index_xy(53, 0)
        peaks.append('Pt111')
        xs[3], ys[3] = data.get_index_xy(39.76, -10)
        xs[4], ys[4] = data.get_index_xy(35, 10)
        xs[5], ys[5] = data.get_index_xy(45, 0)

    elif data.y[0] < 45 < data.y[-1] and data.x[0] < 32 < data.x[-1]:
        peaks.append('110')
        xs[0], ys[0] = data.get_index_xy(32.56, 30)
        xs[1], ys[1] = data.get_index_xy(30, 60)
        xs[2], ys[2] = data.get_index_xy(35, 45)
        peaks.append('111')
        xs[3], ys[3] = data.get_index_xy(40.12, 50)
        xs[4], ys[4] = data.get_index_xy(35, 60)
        xs[5], ys[5] = data.get_index_xy(45, 54.74)
        peaks.append('Pt200_111')
        xs[6], ys[6] = data.get_index_xy(46.24, 45)
        xs[7], ys[7] = data.get_index_xy(41, 65)
        xs[8], ys[8] = data.get_index_xy(51, 54.74)
    """
    peak = '110'
    x1, y1 = data.get_index_xy(32.56, -10)
    x2, y2 = data.get_index_xy(30, 80)
    x3, y3 = data.get_index_xy(35, 45)
    """
    smapT = data.smap.copy().T
    lines = []
    name = os.path.basename(file_list[i])[:23]
    basename = os.path.dirname(file_list[i])
    print(basename)
    print(name)
    for i in range(len(peaks)):
        lines = []
        # get psi lines and fits
        # get y
        lines.append(data.y[ys[i*3]:ys[i*3+1]])
        # psi_int
        lines.append(data.integrate_2d([xs[i*3+1], ys[i*3], xs[i*3+2],
                                        ys[i*3+1]], 'y'))
        # psi_slice
        lines.append(smapT[xs[i*3], ys[i*3]:ys[i*3+1]])

        ret = get_fit(lines[0], lines[1])
        savename = os.path.join(basename, '%s_psi' % peaks[i])
        fits_to_csv2(lines[0], lines[1], name, savename, plot=False)

        # get 2th lines and fits
        # get x
        lines.append(data.x[xs[i*3+1]:xs[i*3+2]])
        # 2th_int
        lines.append(data.integrate_2d([xs[i*3+1], ys[i*3], xs[i*3+2],
                                        ys[i*3+1]], 'x'))
        # line_2th_slice
        lines.append(data.smap[ys[i*3+2], xs[i*3+1]:xs[i*3+2]])

        ret = get_fit(lines[3], lines[4])
        savename = os.path.join(basename, '%s_2th' % peaks[i])
        fits_to_csv2(lines[3], lines[4], name, savename, plot=False)
    # ""


"""Plot the lines
fig = plt.figure()
ax = []
ax.append(fig.add_subplot(121))
ax.append(fig.add_subplot(122))
#ax.append(fig.add_subplot(223))
#ax.append(fig.add_subplot(224))
for i,row in enumerate(psi_int):
    ax[0].plot(lines[0], row, label=names[i])
    #ax[0].plot(lines[0], psi_slice[i], label=names[i])
    ax[1].plot(lines[3], inter_2th[i], label=names[i])
    #ax[1].plot(lines[3], slice_2th[i], label=names[i])
for i, a in enumerate(ax):
    a.set_ylabel('Intensity')
    a.set_xlabel(u'\u03a8[\u00b0]')
    if i<1:
        a.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", ncol=3)
    if i > 0:
        a.set_xlabel(u'2\u03b8[\u00b0]')
plt.tight_layout()
plt.show()
"""


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


if __name__ == '__main__':
    pass
