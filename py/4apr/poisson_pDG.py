from firedrake import *
from firedrake.output import VTKFile

m = 20
mesh = UnitSquareMesh(m, m, quadrilateral=True)

H = FunctionSpace(mesh, 'DG', 0)
print(f'{m} x {m} mesh DG0:')

u = Function(H, name='u')
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
f = Function(H).interpolate(5 * x + 7 * y - 1.11111)
f.rename('f')

x_func = Function(H).interpolate(x)
y_func = Function(H).interpolate(y)
Delta_h = sqrt(jump(x_func)**2 + jump(y_func)**2)

F = jump(u)/Delta_h*jump(v) * dS - f * v * dx

# enforce Dirichlet BCs weakly
h = 1.0 / m
F += (10.0 / h) * u * v * ds((3,4))

solve(F == 0, u, bcs=None,
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

# no obvious way to get boundary flux; following does NOT work
#   n = FacetNormal(mesh)
#   oflux = assemble(dot(grad(u),n) * ds)
uint = assemble(u * dx)
fint = assemble(f * dx)
print(f'  u integral       = {uint:13.6e}')
print(f'  f integral       = {fint:13.6e}')
print('  ... no obvious way to measure conservation success/failure')

VTKFile("result.pvd").write(f,u)
