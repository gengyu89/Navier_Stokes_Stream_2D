"""
CERI 8315 Term Porject: solving continuity and momentum equations for the case
    of constant viscosity - the finite element approach.

Note that pyplot or pylab does not work in the FEniCS environment.
It is suggested that you solve the PDE in one script, deactivate FEniCS and
    plot the solution in another script.
    
Each script must be ran separately and it is not possible to have a main script
    driving the entire project.
"""

from lib.stream_fe import *

# Global variable for solutions' folder.
DIR = 'fe_solution'
print('Output: %s' % DIR)

# Set-up study domain.
M, N = 21, 31  # provide number of vertexes, they will be
xmax, ymax = 1000000, 1500000  # converted into number of elements inside of subroutines

# Default parameters.
eta = 1e+21
g_y = 10      # m/s^2
rho_L = 3200  # left layer
rho_R = 3300  # right layer

# Examine working directory.
outdir = './output/'+DIR+'/pvd/'
csvfdr = './output/'+DIR+'/csv/'
prepare_dir(outdir)
prepare_dir(csvfdr)  # clean up the csv folder
print("")

# Try with default parameters.
filename = 'default'
print("Solving default parameters...")
solve_coupled(M, N, xmax, ymax, rho_L, rho_R, \
              g_y, eta, outdir, filename)
print("")

# Try with different viscosities.
print("Trying different viscosities...")
for eta_i in range(1, 21, 1):  # use a smaller size
	filename = 'viscosity_%02d' % eta_i
	# eta_i *= 1e20
	solve_coupled(M, N, xmax, ymax, rho_L, rho_R, \
	              g_y, eta_i*1e20, outdir, filename)
print("")

# Try with different density constrasts.
print("Trying different densities...")
for constrast in range(100, 1000, 400):
	rho_Li = 3250 - 0.5 * constrast
	rho_Ri = 3250 + 0.5 * constrast
	filename = 'density_%04d' % constrast
	solve_coupled(M, N, xmax, ymax, rho_Li, rho_Ri, \
	              g_y, eta, outdir, filename)
print("")

# Try with different gravity values.
print("Trying different gravities...")
for g_i in range(1, 10, 4):
	filename = 'gravity_%02d' % g_i
	solve_coupled(M, N, xmax, ymax, rho_L, rho_R, \
	              g_i, eta, outdir, filename)
print("")

# Convert the solutions into ASCII format.
pvdpath = outdir
print("Converting file format...")
vtulist = [fname for fname in os.listdir(pvdpath) if fname.endswith('.vtu')]

for fname in vtulist:
	fin = pvdpath + fname
	fout = csvfdr + fname.replace('000000.vtu', '.txt')
	vtk2ascii(fin, fout)

print("Done. De-activate FEniCS and run fe_matplotlib.py for further processing.")
