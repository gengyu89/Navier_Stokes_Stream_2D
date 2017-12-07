"""
CERI 8315 Term Porject: solving velocity field for the case of constant
    viscosity - the Navier-Stokes approach
    
Note that pyplot or pylab does not work in the FEniCS environment.
It is suggested that you solve the PDE in one script, deactivate FEniCS and
    plot the solution in another script.
    
Each script must be ran separately and it is not possible to have a main script
    driving the entire project.
"""

from lib.navier_stokes import *

# Global variable for solutions' folder.
DIR = 'ns_solution'
print('Output: %s' % DIR)

# Scaling factors.
v_mag = 1e-7
p_mag = 1.0

# Default parameters.
eta = 1e+21
g_y = 0.0     # ignore gravity
rho_L = 3200  # left layer
rho_R = 3300  # right layer

# Set-up study domain.
M, N = 21, 31  # provide number of vertexes, convert into number of elements inside of subroutines
rho = 0.5 * (rho_L + rho_R)  # use uniform density for this problem
p, v = -8*p_mag, 0.8*v_mag
xmax, ymax = 1000000, 1500000

# Subfolder for saving relations.
print("Preparing paths for log...")
logpath = './output/'+DIR+'/log/'
# prepare_dir(logpath)  # do not clean up the log
print("")  # it has the plots produced by plot_all_logs.py


# =============================================================================
#  Routine tests
# =============================================================================

# Solve with default parameters.
print("Solving default parameters...")
outdir = './output/'+DIR+'/default/'
prepare_dir(outdir)  # examine output directory
u_max = solve_ns(M, N, xmax, ymax, p, v, \
	rho, g_y, eta, outdir)
print("")  # use an extra line instead of indentation

# Try with different pressures.
# print("Trying different pressures...")
# logfile = logpath + 'pressure.log'
# fp = open(logfile, 'w')

# p_list = {'psc_00': -8, 'psc_03': -8e3, 'psc_06': -8e6, 'psc_09': -8e9}

# for test_name, p_i in p_list.iteritems():
# 	outdir = './output/'+DIR+'/%s/' % test_name
# 	prepare_dir(outdir)  # examine output directory
# 	u_max = solve_ns(M, N, xmax, ymax, p_i, v, \
# 		rho, g_y, eta, outdir)
# 	lout = '%g\t%g\n' % (p_i, u_max)
# 	fp.write(lout)

# fp.close()
# print('File saved as: %s\n' % logfile)

# Try with different velocities.
print("Trying different velocities...")
logfile = logpath + 'velocity.log'
fp = open(logfile, 'w')

for v_i in range(2, 10, 2):
	outdir = './output/'+DIR+'/vel_%02d/' % v_i
	prepare_dir(outdir)  # examine output directory
	u_max = solve_ns(M, N, xmax, ymax, p, v_i/10.0*v_mag, \
		rho, g_y, eta, outdir)
	lout = '%g\t%g\n' % (v_i/10.0*v_mag, u_max)
	fp.write(lout)

fp.close()
print('File saved as: %s\n' % logfile)

print("Done.")
