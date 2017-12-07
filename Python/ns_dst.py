"""
CERI 8315 Term Porject: solving velocity field for the case of constant
    viscosity with Navier-Stokes equation - the double vortex experiment script
    
Nothing special compared with routine tests but there is an additional B.C.
    added at the bottom
"""

from lib.navier_stokes import *

# Global variable for solutions' folder.
DIR = 'ns_tests'
print('Output: %s' % DIR)

# Set-up study domain.
v = 0.8e-7
M, N = 21, 31  # provide number of vertexes, convert into number of elements inside of subroutines
v_b, v_t = 0.8e-7, 0.8e-7
xmax, ymax = 1000000, 1500000  # these affect the Reynold's number

# Default parameters.
eta = 1e+21
g_y = 0.0     # ignore gravity
rho_L = 3200  # left layer
rho_R = 3300  # right layer

# Solve with default parameters.
print("Running default test...")
rho = 0.5 * (rho_L + rho_R)
outdir = './output/'+DIR+'/default/'
prepare_dir(outdir)
double_swirl(M, N, xmax, ymax, v_b, v_t, rho, g_y, eta, outdir)
print("")


# =============================================================================
#  Double swirls tests
# =============================================================================

# Reduce the viscosity to common fluid level.
full_list = {'air': (1.2, 17.4e-6),
             'honey': (1.40e3, 10.0),
             'methanol': (0.792e3, 5.44e-4),
             'glycerol': (1.261e3, 1.2),
             'mercury': (13.534e3, 1.526e-3),
             'granite': (2.70e3, 3e19),
             'dummy_1': (1.153e3, 2.3e13),
             'dummy_2': (1.153e3, 2.3e11),
             'dummy_3': (1.153e3, 2.3e8)}
# structure: {'matter name': (density [kg/m^3], viscosity [Pa.s])}

## TEST I

# Run vortex test with different materials.
sub_list = ['granite', 'dummy_1', 'dummy_2', 'dummy_3']
v_b, v_t = 0.8e-7, 0.8e-7

# Perform single and double swirls experiments.
print("Performing double swirl tests (this will take a while)...")
for matter, data in full_list.iteritems():
    if matter in sub_list:
        rho_i, eta_i = data
        outdir = './output/'+DIR+'/%s/' % matter
        prepare_dir(outdir)
        double_swirl(M, N, xmax, ymax, v_b, v_t, rho_i, g_y, eta_i, outdir)
    else:
        pass

print("")

## TEST II

# Repeat the above process with a longer scale.
sub_list = {'granite': 5.0, 'dummy_1': 50, 'dummy_2': 500, 'dummy_3': 5000}
v_b, v_t = 0.8e-7, 0.8e-7

# Perform single and double swirls experiments.
print("Performing the above tests on a longer scale...")
for matter, data in full_list.iteritems():
    if matter in sub_list:
        rho_i, eta_i = data
        outdir = './output/'+DIR+'/longer_%s/' % matter
        prepare_dir(outdir)
        length = sub_list[matter]
        double_swirl(M, N, xmax, ymax, v_b, v_t, \
            rho_i, g_y, eta_i, outdir, length)
    else:
        pass

print("")
