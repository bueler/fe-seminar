from firedrake import *

mesh = UnitSquareMesh(10,10)

H = FunctionSpace(mesh,'CG',1)
u = Function(H, name='u(x,y)')
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
dsqr = (x - 0.8)**2 + (y - 0.3)**2
f = Function(H).interpolate(exp(-10.0*dsqr))
f.rename('f(x,y)')

F = ( dot(grad(u),grad(v)) - f*v ) * dx

bdry_ids = (1, 2, 3, 4)   # four sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

File("result.pvd").write(f,u)
