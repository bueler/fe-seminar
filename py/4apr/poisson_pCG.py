from firedrake import *
from firedrake.output import VTKFile

m = 20
mesh = UnitSquareMesh(m, m)

k = 1
H = FunctionSpace(mesh, 'CG', k)
print(f'{m} x {m} mesh CG{k}:')

u = Function(H, name='u')
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
f = Function(H).interpolate(5 * x + 7 * y - 1.11111)
f.rename('f')

F = ( dot(grad(u), grad(v)) - f * v ) * dx

# boundary indices on square:                -4-
#   zero Neumann   on sides 1,2            1|   |2
#   zero Dirichlet on top,bottom 3,4         -3-
bdry_ids = (3, 4)
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

# measure conservation failure (much larger than rounding error)
uint = assemble(u * dx)
fint = assemble(f * dx)
n = FacetNormal(mesh)
oflux = assemble(dot(grad(u),n) * ds)
imbalance = - oflux - fint
print(f'  u integral       = {uint:13.6e}')
print(f'  f integral       = {fint:13.6e}')
print(f'  flux out         = {oflux:13.6e}')
print(f'  imbalance        = {imbalance:13.6e}')

VTKFile("result.pvd").write(f,u)
