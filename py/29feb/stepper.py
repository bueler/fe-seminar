# This example solves the time-dependent heat equation
#   u_t = u_xx + u_yy + f(x,y)
# for a mix of Dirichlet and Neumann boundary conditions
# on the unit square.  It demonstrates
#   1. implicit or explicit time-stepping
#   2. conditional stability of explicit time-stepping
#   3. writing time steps to a file
# The problem and the time-stepping are further documented
# in problem.pdf.

from firedrake import *

implicit = True

M = 5
# M = 20
mesh = UnitSquareMesh(M,M)
N, deltat = 50, 0.002
# N, deltat = 5, 0.02

# vector spaces and trial/test functions
H = FunctionSpace(mesh,'CG',1)
unew = Function(H, name='u^n(x,y)')      # initialized to zero
uold = Function(H, name='u^{n-1}(x,y)')  # ditto
v = TestFunction(H)

# set up source functions
x, y = SpatialCoordinate(mesh)
d2 = (x-0.8) * (x-0.8) + (y-0.1) * (y-0.1)
f = Function(H).interpolate(10 * exp(-10 * d2))
f.rename('f(x,y)')
g = Function(H).interpolate(conditional(y>0.5, 2.0, 0.0))
g.rename('g(y) on Gamma_N')

# weak forms differ: where do unew and uold go?
if implicit:
    F = unew * v * dx + deltat * dot(grad(unew), grad(v)) * dx \
        - uold * v * dx
else:
    F = unew * v * dx \
        - uold * v * dx + deltat * dot(grad(uold), grad(v)) * dx
# add sources to weak form
F -= deltat * (f * v * dx + g * v * ds(1))   # side 1 is {x=0} left side

# zero Dirichlet on other three boundaries
bdry_ids = (2, 3, 4)
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

# time-stepping and output
t = 0.0
outfile = File("result.pvd")
outfile.write(unew, time=t)
for j in range(N):
    solve(F == 0, unew, bcs=[BCs],                      # for unew
          solver_parameters = {'snes_type': 'ksponly',
                               'ksp_type': 'preonly',
                               'pc_type': 'lu'})
    uold.assign(unew)
    t += deltat
    outfile.write(unew, time=t)
File("sources.pvd").write(f,g)
