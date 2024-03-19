from firedrake import *
from firedrake.output import VTKFile  # delete this line in older Firedrake versions

lam = 6.0  # parameter in Bratu equation; critical value between 6.0 and 7.0?

M = 100
mesh = UnitSquareMesh(M,M)

H = FunctionSpace(mesh,'CG',1)
u = Function(H, name='u_h(x,y)')
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
d2 = (x - 0.8)**2 + (y - 0.1)**2
f = Function(H).interpolate(exp(-40.0 * d2))   # spot source
f.rename('f(x,y)')

F = ( dot(grad(u),grad(v)) - lam * exp(u) * v - f*v ) * dx

bdry_ids = (1, 2, 3, 4)   # four sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'newtonls',
                           'snes_monitor': None,
                           'snes_converged_reason': None,
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

normu = norm(u, norm_type='L2')
print(f'{M:d} x {M:d} mesh: |u| = {normu:.6f}')

VTKFile("result.pvd").write(f,u)  # older versions of Firedrake use "File()"
