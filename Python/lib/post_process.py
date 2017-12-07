"""
Load the CSV files produced by fe_main.py and perform further processings

This script requires prettytable
    Type
        $ pip install prettytable
    to install it
"""


# =============================================================================
#  Type conversion
# =============================================================================

from __future__ import print_function
import prettytable
import csv
import sys

def csv2table(fin, fout):
    """Easy CSV to ASCII converter."""
    
    table = None
    with open(fin, 'rb') as csvfile:
        content = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in content:
            if table is None:
                table = prettytable.PrettyTable(row)
            else:
                table.add_row(row)
    
    output = open(fout, 'w')
    print(table, file=output)
    
    print("File saved as:", fout)


# =============================================================================
#  Minor subroutines
# =============================================================================

import os
import shutil

def compute_spacing(M, N, xmax, ymax):
    x_interv = float(xmax-0) / (M-1)
    y_interv = float(ymax-0) / (N-1)
    return max(x_interv, y_interv)

def prepare_dir(path):
    """If output directory exists, clean-up previous output;
    otherwise, create a new one."""
    
    if os.path.isdir(path):
        print("Cleaning previous output...")
        shutil.rmtree(path)
        os.mkdir(path)
    else:
        os.mkdir(path)

import numpy as np

def load_ascii(filepath):
    """Easy and efficient way of data reading from plain text,
    with vectorization using Numpy."""
    # filepath <type 'str'>
    #     - please provide full path
    
    dat = np.loadtxt(filepath)
    return dat[:,0], dat[:,1], dat[:,2]

def load_csv(filepath):
    """Read and convert a CSV file in the prettytable format."""
    
    csv = open(filepath, 'r')
    for i in range(2):
        csv.readline()  # skip the first two lines
    
    lin = csv.readline()
    f, p_0, p_1 = [], [], []
    while lin != "":
        if lin.count('+') == 5:
            pass
        else:
            dat = lin.split('|')
            f_i  = float(dat[1].strip())
            p0_i = float(dat[2].strip())
            p1_i = float(dat[3].strip())
            f.append(f_i)
            p_0.append(p0_i)
            p_1.append(p1_i)
        lin = csv.readline()
        
    csv.close()
    return f, p_0, p_1


# =============================================================================
#  Plotting log files
# =============================================================================

from pylab import *
set_printoptions(suppress=True)

import seaborn as sns
sns.set(style="whitegrid")

def fetch_dat(fullpath):
    """Loading a log and converting it into x and y arrays."""
    
    fp = open(fullpath, 'r')
    lin = fp.readline()
    x, y = [], []
    while lin != "":
	data = lin.strip().split()
	x_i, y_i = float(data[0]), float(data[1])
	x.append(x_i)
	y.append(y_i)
	lin = fp.readline()
    fp.close()
    
    # Return in np.array types.
    return array(x), array(y)  # to perform further scaling

def plot_ns(logfile, quantity, unit):
    """Plot and save a figure for the N-S project."""
    
    figpath = logfile.replace('.log', '.png')
    x, u_max = fetch_dat(logfile)
    
    figure
    title('$u_{max}$ vs. %s B.C.' % quantity.lower())
    plot(x, u_max)
    scatter(x, u_max)
    xlabel('%s at the open edge [%s]' % (quantity, unit))
    # ylabel('Maximum velocity magnitude [m/s]')
    show()
    
    savefig(figpath)
    print("File saved as:", figpath)
    close()

def plot_stream(v_mag, logpath_1, logpath_2):
    """A subroutine to plot the solutions from finite-
    element and finite-difference approaches for comparison."""
    
    # v_mag <type 'float'>
    #     - magnifying factor for the velocities
    #     - must be consistent with fd_main.py and fe_matplotlib.py
    # logpath_1 <type 'str'>
    #     - path to the finite-element log file
    # logpath_2 <type 'str'>
    #     - path to the finite-difference log file
    
    # Create plots.
    figure
    lb1, lb2 = 'Finite Element', 'Finite Difference'
    c0, c1 = '#1f77b4', '#ff7f0e'
    
    # Scale all the velocities to 1e-9
    # Scale the viscosity to 1e21
    figpath_1 = logpath_1 + 'fe_vs_fd.png'
    figpath_2 = logpath_2 + 'fe_vs_fd.png'
    
    # Plot density.
    x_1, y_1 = fetch_dat(logpath_1 + 'density.log')
    x_2, y_2 = fetch_dat(logpath_2 + 'density.log')
    
    subplot(2,2,1)
    title('$v_{max}$ vs. density contrast [kg/m$^3$]')
    # scatter(x_1, y_1/1e9, 'k')
    # scatter(x_2, y_2/1e9, 'k')
    plot(x_1, y_1, label=lb1, color=c0)
    plot(x_2, y_2, label=lb2, color=c1)
    # xlabel('[kg/m$^3$]')
    ylabel('%g [m/s]' % v_mag)
    legend()
    
    # Plot gravity.
    x_1, y_1 = fetch_dat(logpath_1 + 'gravity.log')
    x_2, y_2 = fetch_dat(logpath_2 + 'gravity.log')
    
    subplot(2,2,3)
    xlabel('$v_{max}$ vs. gravity [m/s$^2$]')
    # scatter(x_1, y_1/1e9, 'k')
    # scatter(x_2, y_2/1e9, 'k')
    plot(x_1, y_1, label=lb1, color=c0)
    plot(x_2, y_2, label=lb2, color=c1)
    # xlabel('[m/s$^2$]')
    ylabel('%g [m/s]' % v_mag)
    legend()
    
    # Plot viscosity.
    x_1, y_1 = fetch_dat(logpath_1 + 'viscosity.log')
    x_2, y_2 = fetch_dat(logpath_2 + 'viscosity.log')
    
    subplot(2,2,4)
    xlabel('$v_{max}$ vs. viscosity [Pa.s]')
    # scatter(x_1*1e21, y_1/1e9, 'k')
    # scatter(x_2*1e21, y_2/1e9, 'k')
    plot(x_1*1e21, y_1, label=lb1, color=c0)
    plot(x_2*1e21, y_2, label=lb2, color=c1)
    # xlabel('[Pa.s]')
    # ylabel('10$^9$ [m/s]')
    legend()
    
    show()
    # savefig(figpath_1)
    savefig(figpath_2)
    
    # print("File saved as:", figpath_1)
    print("File saved as:", figpath_2)
    close()


# =============================================================================
#  The main subroutine
# =============================================================================

from numba import jit

# Disable the strange warning from contour.py
seterr(divide='ignore', invalid='ignore')

# @jit
def compute_velocity(csvpath, filename, M, N, xmax, ymax, \
    rho_L, rho_R, headline=""):
    """Post-processing for the stream function - finite element approach.
    Basically the same subroutine as the finite difference approach,
    but relies on the file reading from the converted CSV format."""
    
    # path <type 'str'>
    #     - provide incomplete path like './output/fe_sol/csv/'
    #     - for further processing
    
    # Create fig path from csvpath.
    figpath = csvpath.replace('/csv/', '/plots/')
    fullpath = figpath + filename + '.png'
    
    # Load vorticity and stream funtion.
    f_o, p_0, p_1 = \
        load_ascii(csvpath + filename + '_omega.txt')
    f_p, p_0, p_1 = \
        load_ascii(csvpath + filename + '_psi.txt')  # sharing the same p_0, p_1
    
    # Grid the data into matrices.
    x = linspace(0, xmax, M)
    y = linspace(0, ymax, N)
    omega = griddata(p_0, p_1, f_o, x, y, interp='linear')
    psi   = griddata(p_0, p_1, f_p, x, y, interp='linear')
    
    # Create namespace for the finite difference style.
    ynum, xnum = N, M
    ysize, xsize = ymax, xmax
    ystp = float(ysize) / (ynum-1)
    xstp = float(xsize) / (xnum-1)
    
    # Populate dataset.
    rho = zeros([ynum, xnum])
    rho[:,:xnum/2] = rho_L
    rho[:,xnum/2:] = rho_R
    
    # Calculate velocity from stream function
    # following the subroutine stream_2d()
    vx = zeros([ynum, xnum])
    vy = zeros([ynum, xnum])
    for i in range(2, ynum):
        for j in range(2, xnum):
            vx[i-1, j-1] =  (psi[i,j-1] - psi[i-2,j-1]) / 2.0 / ystp
            vy[i-1, j-1] = -(psi[i-1,j] - psi[i-1,j-2]) / 2.0 / xstp
    
    # Compute min and max magnitudes of the velocity field.
    v_mag = sqrt(vx*vx + vy*vy)
    v_max = np.max(v_mag)
    v_min = np.min(v_mag)
    
    # Start creating plots.
    clf()
    fig = figure(1)
    
    suptitle(headline)
    x2plot, y2plot = x/1e3, y/1e3
    colormap = 'viridis_r'
    
    # Plotting vorticity as colormap.
    subplot(2,2,1)
    pc = pcolor(x2plot, y2plot, omega, cmap=colormap)
    title('Vorticity [s$^{-1}$]')
    ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    axis('tight')
    axis('equal')
    cb = colorbar(pc)
    
    # Plotting density as colormap.
    subplot(2,2,2)
    pc = pcolor(x2plot, y2plot, rho, cmap=colormap)  # shading('interp')
    cb = colorbar(pc)
    title('Density [kg/m$^3$]')
    ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    axis('tight')
    axis('equal')
    
    # Plotting stream function as contours.
    subplot(2,2,3)
    CS = contour(x2plot, y2plot, psi, cmap=colormap)
    clabel(CS, inline=1, fontsize=10)
    xlabel('Stream function [km$^2$/s]')
    ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    axis('tight')
    axis('equal')
    cb = colorbar(CS)
    
    # Plotting velocity vector as arrows using internal nodes only.
    subplot(2,2,4)
    qv = quiver(x2plot[1:-1],y2plot[1:-1], \
                vx[1:-1, 1:-1],vy[1:-1, 1:-1], \
                color='C0')  # normalize with average magnitude
    xlabel('Velocity [km/s]')
    ylim([max(y2plot[1:-1]), min(y2plot[1:-1])])  # flip y-axis
    axis('tight')
    axis('equal')
    text(-500, 150, 'Magnitude ($10^{-7}$)')
    text(-500, 300, '[%f, %f]' % (v_min*1e7, v_max*1e7))
    
    show()
    savefig(fullpath)
    print("File saved as:", fullpath)
    
    close()    
    return v_max
