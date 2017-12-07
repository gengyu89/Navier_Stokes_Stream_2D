"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
    - modified for the purpose of term project
Incremental Pressure Correction Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

from __future__ import print_function
from fenics import *
import numpy as np
import time  # newly added


# =============================================================================
#  Minor subroutines
# =============================================================================

import os
import shutil

def prepare_dir(path):
    """If output directory exists, clean-up previous output;
    otherwise, create a new one."""
    if os.path.isdir(path):
        print("Cleaning previous output...")
        shutil.rmtree(path)
        os.mkdir(path)
    else:
        os.mkdir(path)

import sys

def progress_bar(completed, total):
    """The progress bar of FEniCS sucks. Let me use my own."""

    # Check validity.
    assert completed <= total, 'Completed greater than total!'
    
    # Get terminal width.
    terminal_size = os.popen('stty size', 'r').read()
    terminal_width = int(terminal_size.split()[1])
    
    # Remaining and finished bar lengths.
    bar_length = terminal_width - 40
    ratio = float(completed) / float(total)
    finished_length = int(ratio * bar_length)
    remain_length = bar_length - finished_length
    
    # Flashing periods.
    if completed == total:
        nof_period = 0
    else:
        nof_period = completed % 3
    nof_spaces = 3 - nof_period
    cat_period = '.' * nof_period + ' ' * nof_spaces  # the spaces are needed to cover previous screen output
    
    # Prepare status bar and percentage.
    bar = '#' * finished_length + '_' * remain_length
    percentage = ratio * 100
    status = '\r%-2s[%s] %5.1f percent completed.' % (' ', bar, percentage) + cat_period
    
    # Update status.
    sys.stdout.write(status)
    sys.stdout.flush()


# =============================================================================
#  Navier-Stokes and its subroutines
# =============================================================================

def boundary(x, on_boundary):  # not sure if this works
    """Common boundary condition shared by both PDEs."""
    return on_boundary

# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(mu, u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Use numba to accelerate the main subroutine.
# Numba can accelerate all Python codes not limited to numpy!
# from numba import jit, vectorize, float64

# @jit
def solve_ns(M, N, xmax, ymax, p, v, \
    rho, g_y, eta, outdir, length=5.0, slip=False):
    # p <type 'float'>
    #     - pressure at the open edge
    # v <type 'float'>
    #     - right-flow velocity at the open edge
    # length <type 'float'>
    #     - total length in sec, use 50% of the FEniCS example as the default
    # slip <type 'bool'>
    #     False - apply no-slip condition
    #     True  - free-slip
    
    # Using parameters from the example.
    T = float(length)  # final time
    num_steps = 50     # we do not need a very dense time step in this project
    
    if T >= 100:
        num_steps *= 10  # this is to prevent the IPCS from exploding
    else:
        pass
    
    dt = T / num_steps # time step size
    mu = eta           # kinematic viscosity
    K  = 0.13e12       # average bulk modulus in the Earth's upper mantle

    # Estimate Reynold's number using input parameters.
    V = 0.5 * abs(v)       # mean flow velocity
    D = min(xmax, ymax)    # diameter
    Re = rho * V * D / mu  # using the formula for a tube
    print("Reynold's number (estimated): %g" % Re)
    
    # Create mesh and define function spaces
    mesh = RectangleMesh(Point(0, 0), Point(xmax, ymax), M-1, N-1)
    V = VectorFunctionSpace(mesh, 'P', 2)  # velocity has two components
    Q = FunctionSpace(mesh, 'P', 1)
    
    # Define boundaries (this must be modified)
    open_edge = 'near(x[1], %f)' % ymax
    walls = 'near(x[0], 0) || near(x[0], %f) || near(x[1], 0)' % xmax

    # Define boundary conditions
    p_D = Constant(p)  # DirichletBC() does not take variables
    bcp_flow   = DirichletBC(Q, p_D, open_edge)
    u_D = Constant((v, 0))  # write it separately
    bcu_flow   = DirichletBC(V, u_D, open_edge)
    bcu_noslip = DirichletBC(V, Constant((0, 0)), walls)

    bcp = [bcp_flow]
    if slip:
        bcu = [bcu_flow]
    else:
        bcu = [bcu_flow, bcu_noslip]

    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)

    # Define functions for solutions at previous and current time steps
    u_n = Function(V)  # previous
    u_  = Function(V)  # current, u_{n+1}
    p_n = Function(Q)  # previous
    p_  = Function(Q)  # current, p_{n+1}

    # Define expressions used in variational forms
    U   = 0.5*(u_n + u)  # u_{n+1/2} = (u_n + u_{n+1}) / 2
    n   = FacetNormal(mesh)
    f   = Constant((0, -rho*g_y))  # reverse the gravity, paraview cannot flip y-axis
    k   = Constant(dt)
    mu  = Constant(mu)
    rho = Constant(rho)
    
    # print("Advancing time levels...")
    print("Time steps:", num_steps)
    
    # Define variational problem for step 1
    F1 = rho*dot((u - u_n) / k, v)*dx + \
         rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
       + inner(sigma(mu, U, p_n), epsilon(v))*dx \
       + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
       - dot(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)

    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q))*dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

    # Define variational problem for step 3
    a3 = dot(u, v)*dx
    L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

    # Assemble matrices
    A1 = assemble(a1)  # assemble(): used to re-evaluate a matrix expression
    A2 = assemble(a2)  # when you have multiple time-steps so that you do not
    A3 = assemble(a3)  # need to put a long expression into the loop body

    # Apply boundary conditions to matrices
    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]

    # Create VTK file for saving solution
    vtk_u = File(outdir + 'velocity.pvd')
    vtk_p = File(outdir + 'pressure.pvd')

    # Initialize progress bar
    # progress = Progress('Time-stepping')
    # set_log_level(PROGRESS)

    # Time-stepping
    t = 0
    for n in range(num_steps):
        
        # Update current time
        t += dt
        
        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(b1) for bc in bcu]
        solve(A1, u_.vector(), b1)
        
        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcp]
        solve(A2, p_.vector(), b2)
        
        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3)
        
        # Plot solution
        vtk_u << (u_, t)  # updated u
        vtk_p << (p_, t)  # updated p
        # File(outdir + filename + '_rho.pvd') << (rho_, t)  # updated rho
        
        # Compute error
        progress_bar(n+1, num_steps)
        # print('t = %.2f: error = %.3g' % (t, error))
        u_max = u_.vector().array().max()
        # print('max u:', u_max)
        
        # Update previous solution
        u_n.assign(u_)  # current becomes previous
        p_n.assign(p_)  # for the next iteration
        # rho_n = rho_
        
    # Hold plot
    print()
    
    print("Files saved as:", outdir + '*.pvd')
    return u_max  # maximum magnitude from the last time step

# @jit
def double_swirl(M, N, xmax, ymax, vb, v, \
    rho, g_y, eta, outdir, length=5.0, slip=False):
    """A modified version of solve_ns(),
    with an additional B.C. at the bottom,
    two density layers were merged into the same density.
    """
    # vb <type 'float'>
    #     - velocity at the additional edge
    # v  <type 'float'>
    #     - right-flow velocity at the top
    # length <type 'float'>
    #     - total length in sec, use 50% of the FEniCS example as the default
    # slip <type 'bool'>
    #     False - apply no-slip condition
    #     True  - free-slip
    
    # Using parameters from the example.
    T = float(length)  # final time
    num_steps = 50     # we do not need a very dense time step in this project
    
    if T >= 100:
        num_steps *= 10  # this is to prevent the IPCS from exploding
    else:
        pass
    
    dt = T / num_steps # time step size
    mu = eta           # kinematic viscosity
    
    # Estimate Reynold's number using input parameters.
    V = 0.5 * np.mean(abs(vb) + abs(v))  # mean flow speed
    D = min(xmax, ymax)           # diameter
    Re = rho * V * D / mu         # using the formula for a tube
    print("Reynold's number (estimated): %g" % Re)
    
    # Create mesh and define function spaces
    mesh = RectangleMesh(Point(0, 0), Point(xmax, ymax), M-1, N-1)  # totally adjustable
    V = VectorFunctionSpace(mesh, 'P', 2)  # velocity has two components
    Q = FunctionSpace(mesh, 'P', 1)
    
    # Define boundaries (this must be modified)
    bottom, top = 'near(x[1], 0)', 'near(x[1], %f)' % ymax  # you only need to care about
    walls = 'near(x[0], 0) || near(x[0], %f)' % xmax  # your parameters in the arguments
    
    # Define boundary conditions
    bcp_bflow = DirichletBC(Q, Constant(-8), bottom)
    bcp_tflow = DirichletBC(Q, Constant(-8), top)
    u_B = Constant((vb, 0))  # bottom
    u_D = Constant((v, 0))   # top
    bcu_bflow = DirichletBC(V, u_B, bottom)
    bcu_tflow = DirichletBC(V, u_D, top)
    bcu_noslip = DirichletBC(V, Constant((0, 0)), walls)

    bcp = []  # see what happens if we remove the B.C. for pressure
    # bcp = [bcp_bflow, bcp_tflow]
    if slip:
        bcu = [bcu_bflow, bcu_tflow]
    else:
        bcu = [bcu_bflow, bcu_tflow, bcu_noslip]

    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)

    # Define functions for solutions at previous and current time steps
    u_n = Function(V)  # previous
    u_  = Function(V)  # current, u_{n+1}
    p_n = Function(Q)  # previous
    p_  = Function(Q)  # current, p_{n+1}

    # Define expressions used in variational forms
    U   = 0.5*(u_n + u)  # u_{n+1/2} = (u_n + u_{n+1}) / 2
    n   = FacetNormal(mesh)
    f   = Constant((0, -rho*g_y))  # reverse the gravity, paraview cannot flip y-axis
    k   = Constant(dt)
    mu  = Constant(mu)
    rho = Constant(rho)
    
    # print("Advancing time levels...")
    print("Time steps:", num_steps)
    
    # Define variational problem for step 1
    F1 = rho*dot((u - u_n) / k, v)*dx + \
         rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
       + inner(sigma(mu, U, p_n), epsilon(v))*dx \
       + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
       - dot(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)

    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q))*dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

    # Define variational problem for step 3
    a3 = dot(u, v)*dx
    L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

    # Assemble matrices
    A1 = assemble(a1)  # assemble(): used to re-evaluate a matrix expression
    A2 = assemble(a2)  # when you have multiple time-steps so that you do not
    A3 = assemble(a3)  # need to put a long expression into the loop body

    # Apply boundary conditions to matrices
    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]

    # Create VTK file for saving solution
    vtk_u = File(outdir + 'velocity.pvd')
    vtk_p = File(outdir + 'pressure.pvd')

    # Time-stepping
    t = 0
    for n in range(num_steps):
        
        # Update current time
        t += dt
        
        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(b1) for bc in bcu]
        solve(A1, u_.vector(), b1)
        
        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcp]
        solve(A2, p_.vector(), b2)
        
        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3)
        
        # Plot solution
        vtk_u << (u_, t)  # updated u
        vtk_p << (p_, t)  # updated p
        # File(outdir + filename + '_rho.pvd') << (rho_, t)  # updated rho
        
        # Compute error
        progress_bar(n+1, num_steps)
        # print('t = %.2f: error = %.3g' % (t, error))
        u_max = u_.vector().array().max()
        # print('max u:', u_max)
        
        # Update previous solution
        u_n.assign(u_)  # current becomes previous
        p_n.assign(p_)  # for the next iteration
        # rho_n = rho_
        
    # Hold plot
    print()
    
    print("Files saved as:", outdir + '*.pvd')
    return u_max  # maximum magnitude from the last time step
