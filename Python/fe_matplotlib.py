"""
The stream function - vorticity approach
    - script for further processing: compute velocities and plot

Please run in Wing IDE
    - matplotlib has problem with figure saving and closing
"""

from lib.post_process import *

# Global variable for solutions' folder.
DIR = 'fe_solution'
print('Output: %s' % DIR)

# Set-up study domain.
M, N = 21, 31  # for Python, provide the number of grids
xmax, ymax = 1000000, 1500000

# Default parameters.
eta = 1e+21
g_y = 10      # m/s^2
rho_L = 3200  # left layer
rho_R = 3300  # right layer

# Examine working directories.
logpath = './output/'+DIR+'/log/'  # output
figpath = './output/'+DIR+'/plots/'  # output
prepare_dir(logpath)  # clean up the log
prepare_dir(figpath)
print("")

# Convert CSV into TXT files.
csvpath = './output/'+DIR+'/csv/'  # input
csvlist = [fname for fname in os.listdir(csvpath) \
           if fname.endswith('.csv')]

# print("Converting file format...")
# for csv in csvlist:
#     fin = csvpath + csv  # full path with type extension
#     fout = fin.replace('.csv', '.txt')
#     csv2table(fin, fout)  # call csv2table() to perform conversion
# print("")


# =============================================================================
#  Restart the routine processes
# =============================================================================

import time
v_mag = 1e7  # magnifying factor for velocities in the log files

## The following processes are basically the same as the main process.
## It is just that Matplotlib has to work in a separate script.

# Process default parameters.
filename = 'default'
print("Processing default parameters...")
compute_velocity(csvpath, filename, M, N, \
                xmax, ymax, rho_L, rho_R)
print("")

# Process different viscosities.
print("Processing different viscosities...")
logfile = logpath + 'viscosity.log'
fp = open(logfile, 'w')

print("<Performance test>")
start = time.time()

for eta_i in range(1, 21, 1):  # must be consistent with fe_main.py
    filename = 'viscosity_%02d' % eta_i
    headline = 'Viscosity (10$^{21}$): %f Pa.s' % (eta_i/10.0)
    u_max = compute_velocity(csvpath, filename, M, N, \
                            xmax, ymax, rho_L, rho_R, headline)
    lout = '%f\t%f\n' % (eta_i/10.0, u_max*v_mag)
    fp.write(lout)
    
end = time.time()
print('<Elapsed time> %f sec' % (end-start))
    
fp.close()
print("File saved as: %s\n" % logfile)

# Process different density constrasts.
print("Processing different densities...")
logfile = logpath + 'density.log'
fp = open(logfile, 'w')

for constrast in range(100, 1000, 400):
    rho_Li = 3250 - 0.5 * constrast
    rho_Ri = 3250 + 0.5 * constrast
    headline = 'Density contrast: %d kg/m$^3$' % constrast
    filename = 'density_%04d' % constrast
    u_max = compute_velocity(csvpath, filename, M, N, \
                            xmax, ymax, rho_Li, rho_Ri, headline)
    lout = '%f\t%f\n' % (constrast, u_max*v_mag)
    fp.write(lout)
    
fp.close()
print("File saved as: %s\n" % logfile)

# Process different gravity values.
print("Processing different gravities...")
logfile = logpath + 'gravity.log'
fp = open(logfile, 'w')

for g_i in range(1, 10, 4):
    headline = '$g$ = %f km/s$^2$' % g_i
    filename = 'gravity_%02d' % g_i
    u_max = compute_velocity(csvpath, filename, M, N, \
                            xmax, ymax, rho_L, rho_R, headline)
    lout = '%f\t%f\n' % (g_i, u_max*v_mag)
    fp.write(lout)
    
fp.close()
print("File saved as: %s\n" % logfile)

import warnings
# warnings.simplefilter('ignore')

print 'Done. Velocities are in the scale of %g in the log files.' % v_mag
