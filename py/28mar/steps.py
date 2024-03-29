from firedrake import *
from firedrake.petsc import PETSc
from firedrake.output import VTKFile  # in older Firedrake versions: delete this line
import numpy as np

# Setting up the Mesh.
m = 50   # regular m x m mesh

mesh = UnitSquareMesh(m, m)

# Get mesh information
mesh.topology_dm.viewFromOptions('-dm_view')

# We'll define our obstacle and solution space to be CG1
V = FunctionSpace(mesh, "CG", 1)

# Load will be pushing down on the surface of the membrane
f = Constant(-20)

u = Function(V, name="Solution")
v = TestFunction(V)


# Objective for poisson equation with source term f
J = 0.5*inner(grad(u), grad(u))*dx - inner(f, u)*dx
# Derivative of the objective, gives the VI residual. F >= 0
F = derivative(J, u, v)

# Zero boundary conditions on the boundaries labeled 1 and 2.
bcs = DirichletBC(V, 0, (1, 2))


# Define the obstacle function psi
x = SpatialCoordinate(mesh)[0]
obstacle = conditional(lt(x, +0.25), -0.2,
                       conditional(lt(x, +0.50), -0.4,
                       conditional(lt(x, +0.75), -0.6,
                                   -0.8)))

lb = Function(V).interpolate(obstacle)
lb.rename("LowerBound")


# All the same solver parameters as before.

sp = {"snes_vi_monitor": None,         # prints residual norms and # of active set nodes for each Newton iteration
      "snes_type": "vinewtonrsls",  # VI Newton with with reduced space line search
      "snes_converged_reason": None,  # prints CONVERGED_... message at end of solve
      "snes_rtol": 1.0e-8,  # Relative convergence tolerance
      # Absolute convergence tolerance (Applied to the residual of basis functions on the inactive set)
      "snes_atol": 1.0e-12,
      "snes_stol": 1.0e-12,  # Norm of step convergence tolerance
      "snes_vi_zero_tolerance": 1.0e-12,
      "snes_linesearch_type": "basic",
      # these 3 options say Newton step equations are solved by LU
      "ksp_type": "preonly",
      "pc_type": "lu",
      "pc_factor_mat_solver_type": "mumps"}
problem = NonlinearVariationalProblem(F, u, bcs)
solver = NonlinearVariationalSolver(
    problem, solver_parameters=sp, options_prefix="")
ub = Function(V).interpolate(Constant(PETSc.INFINITY))  # no upper obstacle
solver.solve(bounds=(lb, ub))


outfile = 'Steps.pvd'
if outfile:
    PETSc.Sys.Print('writing to %s ...' % outfile)
    towrite = (u,)
    lb.rename("lb (lower obstacle)")
    lgap = Function(V).interpolate(u - lb)
    lgap.rename("lgap = u - lb")
    towrite += (lb, lgap)
    # or File(..) in older Firedrake:
    VTKFile(outfile).write(*towrite)  # unpack tuple as arguments
