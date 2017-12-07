from __future__ import print_function  # required by FEniCS
from fenics import *


# =============================================================================
#  Minor subroutines
# =============================================================================

import os
import shutil

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """Float number comparison with tolerance, shared by both PDEs."""
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def boundary(x, on_boundary):  # not sure if this works
    """Common boundary condition shared by both PDEs."""
    return on_boundary

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
#  Type conversion
# =============================================================================

import numpy as np
np.set_printoptions(suppress=True)

def split_data(lin):
    """Since the data size is huge,
    use numpy to accelerate it."""
    
    split_1 = lin.split('<')
    str_1 = split_1[1]
    split_2 = str_1.split('>')
    str_2 = split_2[1]
    
    return np.array(str_2.split())  # return a list of splitted strings

def ispoints(lin):
    dat = lin.strip()
    key_1 = 'type="Float64"'
    key_2 = 'NumberOfComponents="3"'
    key_3 = 'format="ascii"'
    return dat.endswith('</DataArray>') and \
           key_1 in dat and key_2 in dat and key_3 in dat

def isvalue(lin):
    dat = lin.strip()
    key_1 = 'type="Float64"'
    key_2 = 'format="ascii"'
    return dat.endswith('</DataArray>') and \
           key_1 in dat and key_2 in dat

def vtk2ascii(vtk, ascii):
    """Data size is huge. Properly use
    vectorization to accelerate the processes."""
    # vtk, ascii <type 'str'> - full paths
    
    # Read data.
    fin = open(vtk, 'r')
    lin = fin.readline()
    c_0, c_1 = 0, 0  # set two counters
    
    while lin != "":
        if ispoints(lin):
            c_0 += 1  # after the loop, c_0 must be exactly one
            p = split_data(lin)
            size = len(p)
            assert size % 3 == 0, \
                   "Not divisible by three! Invalid data size."
            p_reshaped = p.reshape(size/3, 3)
            p_0, p_1 = p_reshaped[:,0], p_reshaped[:,1]  # np.array
        elif isvalue(lin):
            c_1 += 1  # after the loop, c_1 must be exactly one
            f = split_data(lin)  # np.array
        else:
            pass
        lin = fin.readline()
        
    fin.close()
    
    # Check validity.
    assert c_0 == 1, \
           "Points reading duplicated! Invalid data processing."
    assert c_1 == 1, \
           "Values reading duplicated! Invalid data processing."
    assert len(p_0) == len(f), \
           "Vector lengths do not match! Data reading failed."
    
    # Write data into file.
    fout = open(ascii, 'w')
    
    for i, f_i in enumerate(f):
        p_0i, p_1i = p_0[i], p_1[i]
        lout = '%s\t%s\t%s\n' % (f_i, p_0i, p_1i)  # kept the original 'str' type
        fout.write(lout)
        
    fout.close()
    print("File saved as:", ascii)


# =============================================================================
#  Solve the PDE system
# =============================================================================

def solve_vorticity(M, N, xmax, ymax, rho_L, rho_R, g_y, eta, outdir, filename):
    """Solve the first equation"""
    
    # M, N, xmax, ymax - essential parameters to create mesh
    # rho_L, rho_R, g_y, eta - parameters for calculating force function
    
    # Set-up study domain.
    lin_x = np.linspace(0, xmax, M)
    mid_1, mid_2 = lin_x[M/2-1], lin_x[M/2+1]
    # print("Midpoints:", mid_1, mid_2)
    x_spacing = float(xmax) / (M-1)  # convert vertexes into number of intervals
    
    # Create mesh and define function space
    mesh = RectangleMesh(Point(0.0, 0.0), Point(xmax, ymax), M-1, N-1)
    V = FunctionSpace(mesh, 'P', 1)  # convert vertexes into number of elements
    
    # Define boundary condition
    u_D = Constant(0.0)  # what does degree mean?
    bc = DirichletBC(V, u_D, boundary)  # call boundary()
    
    # Define force function
    peak = g_y/eta * (rho_R-rho_L)/(2.0*x_spacing)  # use central difference representation
    force = Expression('mid_1<=x[0] && x[0]<=mid_2 ? peak : 0', \
        degree=2, mid_1=mid_1, mid_2=mid_2, peak=peak)
    
    # use consistent degree as the element type in general, according to the tutorial
    # use a higher degree if it is an exact solution for computing error
    
    # Define variational problem
    omega = TrialFunction(V)
    v = TestFunction(V)
    f = force
    a = -dot(grad(omega), grad(v))*dx  # even in 2-D, this is still dx
    L = f*v*dx
    
    # Compute solution
    omega = Function(V)
    solve(a == L, omega, bc)
    
    # Save intermediate files.
    File(outdir + filename + '_mesh.xml') << mesh  # XML for intermediate results
    File(outdir + filename + '_omega.xml') << omega  # VTK for final solutions
    
    # Save solution to file in VTK format
    # File(outdir + filename + '_rho.pvd') << rho  # only solutions can be written into files
    File(outdir + filename + '_omega.pvd') << omega
    
    # Hold plot
    interactive()
    
def solve_stream(outdir, filename):
    """Solve the second equation"""
    
    # Create mesh and define function space
    mesh = Mesh(outdir + filename + '_mesh.xml')  # load mesh from previous session
    V = FunctionSpace(mesh, 'P', 1)  # must be consistent with vorticity
    
    # Define boundary condition
    u_D = Constant(0.0)
    bc = DirichletBC(V, u_D, boundary)
    
    # Define variational problem
    psi = TrialFunction(V)
    v = TestFunction(V)
    f = Function(V, outdir + filename + '_omega.xml')
    a = -dot(grad(psi), grad(v))*dx
    L = f*v*dx
    
    # Compute solution
    psi = Function(V)
    solve(a == L, psi, bc)
    
    # Save solution to file in VTK format
    File(outdir + filename + '_psi.pvd') << psi
    
    # Throw intermediate files.
    os.remove(outdir + filename + '_omega.xml')
    os.remove(outdir + filename + '_mesh.xml')
    
    # Hold plot
    interactive()
    
def solve_coupled(M, N, xmax, ymax, rho_L, rho_R, g_y, eta, outdir, filename):
    """Solve the PDE system.
    Driver script for the above two subroutines."""
    
    solve_vorticity(M, N, xmax, ymax, rho_L, rho_R, g_y, eta, outdir, filename)
    solve_stream(outdir, filename)

    print("Files saved as:", outdir+filename+'_*.pvd')
