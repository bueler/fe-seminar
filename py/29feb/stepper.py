# This example shows:
#   1. implicit or explicit time-stepping
#   2. conditional stability of explicit time-stepping
#   3. writing time steps to a file
# The PDE problem and the time-stepping are documented
# in problem.pdf.

from firedrake import *

M = 5
# M = 50
mesh = UnitSquareMesh(M,M)
N, deltat = 50, 0.001
# N, deltat = 5, 0.01
implicit = True

H = FunctionSpace(mesh,'CG',1)
unew = Function(H, name='u^n(x,y)')      # initialized to zero
uold = Function(H, name='u^{n-1}(x,y)')  # ditto
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
f = Function(H).interpolate(Constant(0.0)) # or something else ...
f.rename('f(x,y)')
g = Function(H).interpolate(Constant(1.0))
g.rename('g on Gamma_N')

if implicit:
    F = (unew * v + deltat * dot(grad(unew), grad(v))  ) * dx \
        - ((uold + deltat * f) * v) * dx \
        - deltat * g * v * ds(1)
else:
    F = (unew * v) * dx \
        - ((uold + deltat * f) * v - deltat * dot(grad(uold), grad(v))) * dx \
        - deltat * g * v * ds(1)

bdry_ids = (2, 3, 4)   # other three sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

# time-stepping loop
t = 0.0
outfile = File("result.pvd")
outfile.write(unew, time=t)
for j in range(N):
    solve(F == 0, unew, bcs=[BCs],
          solver_parameters = {'snes_type': 'ksponly',
                               'ksp_type': 'preonly',
                               'pc_type': 'lu'})
    uold.assign(unew)
    t += deltat
    outfile.write(unew, time=t)
