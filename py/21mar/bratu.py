# This example solves the nonlinear and steady-state Liouville-Bratu equation
#   - u_xx - u_yy - lambda e^u = 0
# on the unit square, with zero Dirichlet boundary conditions.  The equation
# is a model for thermal runaway, for instance in chemical explosives; see
#   https://en.wikipedia.org/wiki/Liouville%E2%80%93Bratu%E2%80%93Gelfand_equation
# This example illustrates that nonlinear PDEs can have critical parameter values
# where solutions cease to exist. This code demonstrates basic Newton iteration,
# and demonstrates some corresponding PETSc SNES options.
# Documentation:  slides/21mar-bratu.pdf

from firedrake import *
from firedrake.output import VTKFile  # in older Firedrake versions: delete this line

lam = 6.0   # parameter in Bratu equation; critical value occurs between 6.0 and 7.0

M = 20      # note that the precise critical value is affected by resolution
            # ... consider larger M if you want to pin down the critical value

mesh = UnitSquareMesh(M,M)
H = FunctionSpace(mesh,'CG',1)
u = Function(H, name='u_h(x,y)')
v = TestFunction(H)

# nonlinear weak form with "exp(u)"
F = ( dot(grad(u),grad(v)) - lam * exp(u) * v ) * dx

# optionally add a source term, for example:
# x, y = SpatialCoordinate(mesh)
# d2 = (x - 0.8)**2 + (y - 0.1)**2
# f = Function(H).interpolate(exp(-40.0 * d2))   # spot source
# f.rename('f(x,y)')
# F -= f * v * dx

bdry_ids = (1, 2, 3, 4)   # four sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

solve(F == 0, u, bcs=[BCs],
      options_prefix='s', # add -s_xxx_yyy options at runtime
      solver_parameters = {'snes_type': 'newtonls',
                           'snes_linesearch_type': 'basic', # or: 'bt'
                           'snes_monitor': None,
                           'snes_converged_reason': None,
                           'ksp_type': 'preonly', # or cg, etc.
                           'pc_type': 'lu'}) # or multigrid, etc.

normu = norm(u, norm_type='L2')
print(f'{M:d} x {M:d} mesh: |u| = {normu:.6f}')

print('writing solution to result.pvd ...')
VTKFile("result.pvd").write(u)  # older versions of Firedrake use "File()"
