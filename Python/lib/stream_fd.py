from __future__ import print_function


# =============================================================================
#  Basic subroutines
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


# =============================================================================
#  Mesh grid and force function
# =============================================================================

from pylab import *
set_printoptions(suppress=True)

import seaborn as sns
sns.set(style="whitegrid")

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

def plot_meshgrids(M, N, xmax, ymax, outdir):
    """A function to show mesh grids of the study domain,
    used for quadrilateral elements only."""
    
    # Compute spacing.
    h = compute_spacing(M, N, xmax, ymax)
    
    # Set-up study domain.
    x = linspace(0, xmax, M)  # in [m]
    y = linspace(0, ymax, N)  # in [m]
    
    figure
    title('Mesh Grid')
    
    # Plot horizontal lines.
    for y_i in y:
        x_to_plot = array([0, xmax])
        y_to_plot = array([y_i, y_i])
        plot(x_to_plot/1e3, y_to_plot/1e3)
    
    # Plot vertical lines.
    for x_i in x:
        x_to_plot = array([x_i, x_i])
        y_to_plot = array([0, ymax])
        plot(x_to_plot/1e3, y_to_plot/1e3)
    
    # Flip y-axis.
    plt.gca().invert_yaxis()
    
    xlabel('$x$ [km]')
    ylabel('$y$ [km]')
    axis('equal')
    text(-500, 0, 'h = %f km' % (h/1e3))
    
    show()
    fullpath = outdir + 'mesh_grid.png'
    savefig(fullpath)
    print('File saved as: %s' % fullpath)
    close()

def force_function(M, N, xmax, ymax, rho_L, rho_R, g_y, eta, outdir):
    
    fullpath = outdir + 'force_function.png'
    
    # Set-up study domain.
    x = linspace(0, xmax, M)  # in [m]
    y = linspace(0, ymax, N)  # in [m]
    
    # Populate dataset.
    rho = zeros([N, M])
    rho[:,:M/2] = rho_L  # M/2 is an integer
    rho[:,M/2:] = rho_R
    
    f = zeros([N, M])  # the force function
    f[:,0] = 0.0  # this can be actually ignored
    f[:,-1] = 0.0
    
    for j in range(1, M-1):
        f[:,j] = (rho[:,j+1]-rho[:,j-1]) / (x[j+1]-x[j-1])
    
    f_xy = f * g_y / eta  # normalize the entire matrix
    
    # Visualize the force function
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.gca(projection='3d')  # make sure Axes3D is imported
    title('Force Function')
    
    # Grid the data.
    x_grid, y_grid = np.meshgrid(x, y)
    sfc = ax.plot_surface(x_grid/1e3, y_grid/1e3, f_xy, \
                          cmap='viridis_r', linewidth=0.0, antialiased=False)
    
    # Flip y-axis.
    plt.gca().invert_yaxis()
    plt.colorbar(sfc)
    
    ax.view_init(30, -30)  # default: (30, -60) x closer
        # we make y closer so that it looks longer than x
    xlabel('$x$ [km]')
    ylabel('$y$ [km]')
    
    show()
    savefig(fullpath)
    print('File saved as: %s' % fullpath)
    close()


# =============================================================================
#  The stream function
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# Block the stupid warning from contour.py
np.seterr(divide='ignore', invalid='ignore')

from numba import jit, vectorize, float64

# @jit #(nopython=True)
def stream_2d(xsize, ysize, xnum, ynum, eta, gy, rho_L, rho_R, \
              outdir, filename='default', headline=''):
    """Solution of 2D Stokes and continuity equations with finite differences
    on a regular grid using stream function - vorticity formulation
    for a medium with constant viscosity"""
    
    # fullpath = '%s%s.png' % (outdir, filename)  # numba does not allow string
    fullpath = '{0}{1}.png'.format(outdir, filename)  # concatenation using '+'
    
    # Grid step.
    xstp = float(xsize) / (xnum-1)  # horizontal
    ystp = float(ysize) / (ynum-1)  # vertical
    
    # Making vectors for nodal points positions.
    x = linspace(0, xsize, xnum)  # Horizontal
    y = linspace(0, ysize, ynum)  # Vertical
    
    # Populate dataset.
    rho = zeros([ynum, xnum])
    rho[:,:xnum/2] = rho_L  # xnum/2 is an integer
    rho[:,xnum/2:] = rho_R
    
    # Matrix of coefficients initialization.
    L = zeros([ynum*xnum, ynum*xnum])
    
    # Vector of right part initialization.
    R = zeros([ynum*xnum, 1])
    
    # Solving Poisson equation for vorticity.
    for i in range(1, ynum+1):
        for j in range(1, xnum+1):
            # Global index for current node.
            k = (j-1) * ynum + i
            # Boundary nodes.
            if i==1 or i==ynum or j==1 or j==xnum:
                # print "if executed!"
                L[k-1,k-1] = 1.0
                R[k-1,0] = 0.0
            # Internal nodes.
            else:
                # Left part.
                L[k-1, k-ynum-1] = 1.0 / (xstp*xstp)
                L[k-1, k-2] = 1.0 / (ystp*ystp)
                L[k-1, k-1] = -2.0 / (xstp*xstp) - 2.0 / (ystp*ystp)
                L[k-1, k] = 1.0 / (ystp*ystp)
                L[k-1, k+ynum-1] = 1.0 / (xstp*xstp)
                # Right part.
                R[k-1, 0] = gy/eta * (rho[i-1,j]-rho[i-1,j-2])/2.0/xstp
    
    # Obtaining vector of solutions S[]
    S = inv(L) * np.matrix(R)
    # print S*1e12
    
    # Reload solutions S[] to 2D vorticity array.
    OMEGA = zeros([ynum, xnum])
    # Process all Grid points.
    for i in range(1, ynum+1):
        for j in range(1, xnum+1):
            # Global index for current node.
            k = (j-1) * ynum + i
            OMEGA[i-1, j-1] = S[k-1]  # no idea why they swap it here
    
    # print OMEGA*1e12
    
    # Solving Poisson equation for stream function.
    R = S  # creates right parts from previous solution
    
    # Obtaining vector of solutions S[]
    S = inv(L) * R
    # print S
    
    # Reload solutions S[] to 2D stream function array PSI[]
    PSI = zeros([ynum, xnum])
    # Process all grid points.
    for i in range(1, ynum+1):
        for j in range(1, xnum+1):
            # Global index for current node.
            k = (j-1) * ynum + i
            PSI[i-1, j-1] = S[k-1]  # no idea why they swap it here
    
    # print PSI
    
    # Compute vx, vy for internal nodes.
    vx = zeros([ynum, xnum])
    vy = zeros([ynum, xnum])
    # Process internal Grid points.
    for i in range(2, ynum):
        for j in range(2, xnum):
            vx[i-1, j-1] =  (PSI[i,j-1] - PSI[i-2,j-1]) / 2.0 / ystp
            vy[i-1, j-1] = -(PSI[i-1,j] - PSI[i-1,j-2]) / 2.0 / xstp
    
    # Compute min and max magnitudes of the velocity field.
    v_mag = sqrt(vx*vx + vy*vy)
    # print v_mag * 1e9
    v_max = np.max(v_mag)
    v_min = np.min(v_mag)
    v_mean = mean(v_mag)
    
    # Plotting solution.
    # Making new figure.
    plt.clf()  # clean-up previous output
    fig = plt.figure(1)
    plt.suptitle(headline)
    x2plot, y2plot = x/1e3, y/1e3  # convert into [km]
    
    colormap = 'viridis_r'  # hot material going up, colormap must be monotonic
    # other possible options: jet, magma, inferno, Spectral, coolwarm
    # do not use: rainbow, hsv, gist_rainbow, etc.
    
    # Plotting vorticity as colormap.
    plt.subplot(2,2,1)
    pc = plt.pcolor(x2plot, y2plot, OMEGA, cmap=colormap)
    plt.title('Vorticity [s$^{-1}$]')
    plt.ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    plt.axis('tight')
    plt.axis('equal')
    cb = plt.colorbar(pc)
    # cb.ax.set_title('[s$^{-1}$]')
    
    # Plotting density as colormap.
    plt.subplot(2,2,2)
    pc = plt.pcolor(x2plot, y2plot, rho, cmap=colormap)  # shading('interp')
    cb = plt.colorbar(pc)
    # cb.ax.set_title('[kg/m$^3$]')  # cb.set_label('[kg/m^3]', rotation=270)
    plt.title('Density [kg/m$^3$]')
    plt.ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    plt.axis('tight')
    plt.axis('equal')
    
    # Plotting stream function as contours.
    plt.subplot(2,2,3)
    CS = plt.contour(x2plot, y2plot, PSI, cmap=colormap)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('Stream function [km$^2$/s]')
    plt.ylim([max(y2plot), min(y2plot)])  # reverse y-axis
    plt.axis('tight')
    plt.axis('equal')
    cb = plt.colorbar(CS)
    # cb.ax.set_title('[km$^2$/s]')
    
    # Plotting velocity vector as arrows using internal nodes only.
    plt.subplot(2,2,4)
    qv = plt.quiver(x2plot[1:-1],y2plot[1:-1], \
                vx[1:-1, 1:-1],vy[1:-1, 1:-1], \
                color='C0')  # normalize with average magnitude
    plt.xlabel('Velocity [km/s]')
    plt.ylim([max(y2plot[1:-1]), min(y2plot[1:-1])])  # reverse y-axis
    plt.axis('tight')
    plt.axis('equal')
    plt.text(-500, 150, 'Magnitude ($10^{-7}$)')
    plt.text(-500, 300, '[%f, %f]' % (v_min*1e7, v_max*1e7))
    
    plt.show()
    plt.savefig(fullpath)
    print('File saved as: %s' % fullpath)  # jit requires Python 3 style
    
    plt.close()
    return v_max