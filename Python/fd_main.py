"""
Main script of the finite difference approach
    with acceleration using Numba

Created on: 2017-11-20
    GENG, Yu
    
Do not run in Terminal
    - the Matplotlib window has a problem with closing
"""

import os
import shutil

# from show_mesh import *
from lib.stream_fd import *

# Global variable for solutions' folder.
DIR = 'fd_solution'
print('Output: %s' % DIR)

# Set-up study domain.
M, N = 21, 31  # for Python, provide the number of grids
xmax, ymax = 1000000, 1500000

# Default parameters.
eta = 1e+21
gy  = 10      # m/s^2
rho_L = 3200  # left layer
rho_R = 3300  # right layer

# Examine working directory.
figpath = './output/'+DIR+'/plots/'
logpath = './output/'+DIR+'/log/'
prepare_dir(figpath)
# prepare_dir(logpath)  # do not clean up the log
print("")  # it has the plots produced by plot_all_logs.py

# Visualize the grid and the force function.
print "Computing grid mesh and force function..."
plot_meshgrids(M, N, xmax, ymax, figpath)  # using default parameters
force_function(M, N, xmax, ymax, rho_L, \
               rho_R, gy, eta, figpath)
print("")


# =============================================================================
#  Start solving the equations
# =============================================================================

import time
v_mag = 1e7  # magnifying factor for velocities in the log files

# Plot with default parameters.
print("Trying default parameters...")
v_max = stream_2d(xmax, ymax, M, N, eta, gy, \
                  rho_L, rho_R, figpath)
print("")

# Try with different viscosities.
print("Trying different viscosities...")
logfile = logpath + 'viscosity.log'
fp = open(logfile, 'w')

# print("<Performance Test>")
# start = time.time()

for eta_i in range(1, 21, 1):  # use a smaller stepsize to track the behavior of viscosity
    if eta_i == 2:
        print("<Performance Test>")
        start = time.time()
    else:
        pass
    filename = 'viscosity_%02d' % eta_i
    suptitle = 'Viscosity (10$^{21}$): %f Pa.s' % (eta_i/10.0)
    # eta_i *= 1e20  # it will break the traverser
    # do not try to mutate the variable being updated in a for loop
    v_max = stream_2d(xmax, ymax, M, N, eta_i*1e20, gy, \
                      rho_L, rho_R, figpath, filename, suptitle)
    lout = '%f\t%f\n' % (eta_i/10.0, v_max*v_mag)    
    fp.write(lout)

end = time.time()
print('<Elapsed Time> %f sec' % (end-start))

fp.close()
print('File saved as: %s\n' % logfile)

# Try with different density constrasts.
print("Trying different densities...")
logfile = logpath + 'density.log'
fp = open(logfile, 'w')

for constrast in range(100, 1000, 200):
    rho_Li = 3250 - 0.5 * constrast
    rho_Ri = 3250 + 0.5 * constrast
    suptitle = 'Density contrast: %d kg/m$^3$' % constrast
    filename = 'density_%04d' % constrast
    v_max = stream_2d(xmax, ymax, M, N, eta, gy, \
                      rho_Li, rho_Ri, figpath, filename, suptitle)
    lout = '%f\t%f\n' % (constrast, v_max*v_mag)    
    fp.write(lout)
    
fp.close()
print('File saved as: %s\n' % logfile)

# Try with different gravity values.
print("Trying different gravities...")
logfile = logpath + 'gravity.log'
fp = open(logfile, 'w')

for g_i in range(1, 10, 2):
    suptitle = '$g$ = %f km/s$^2$' % g_i
    filename = 'gravity_%02d' % g_i
    v_max = stream_2d(xmax, ymax, M, N, eta, g_i, \
                      rho_L, rho_R, figpath, filename, suptitle)
    lout = '%f\t%f\n' % (g_i, v_max*v_mag)
    fp.write(lout)
    
fp.close()
print('File saved as: %s\n' % logfile)

import warnings
print('Done. Velocities are in the scale of %g in the log files.' % v_mag)
