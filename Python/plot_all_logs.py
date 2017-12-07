"""
Script for visualizing all the log files from the Navier-Stokes and the Stream
    function - finite element approach projects

Everything created by this script should be saved in a higher-level directory
    than the raw data, since re-running the computation scripts will clean-up
    the data folder.

Please run this script with Wing IDE. Matplotlib has troubles with saving and
    closing figures in the Terminal.
"""

from lib.post_process import *

# Set paths to dependencies.
NSDIR = 'ns_solution'
FEDIR = 'fe_solution'
FDDIR = 'fd_solution'


# =============================================================================
#  Plot the Navier-Stokes project
# =============================================================================

# Path to the log files
# must be consistent with ns_main.py
logpath = './output/'+NSDIR+'/log/'

# Plot the pressure log.
# logfile = logpath + 'pressure.log'
# quantity, unit = 'Pressure', 'Pa'
# plot_ns(logfile, quantity, unit)

# Plot the velocity log.
logfile = logpath + 'velocity.log'
quantity, unit = 'Velocity', 'm/s'
plot_ns(logfile, quantity, unit)


# =============================================================================
#  Plot the stream function - finite element approach
# =============================================================================

# Paths to the log files
# must be consistent with fe_matplotlib.py and fd_main.py
logpath_1 = './output/'+FEDIR+'/log/'
logpath_2 = './output/'+FDDIR+'/log/'

# Call subroutine to create the plots.
v_mag = 1e7  # must be consistent with fd_main.py and fe_matplotlib.py
plot_stream(v_mag, logpath_1, logpath_2)
print("Done.")
