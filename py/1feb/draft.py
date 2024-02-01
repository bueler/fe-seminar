from firedrake import *

# TODO:  put loop around the solver, with increasing resolution,
#        and measure convergence
mesh = UnitSquareMesh(5,5)

H = FunctionSpace(mesh,'CG',1)
u = Function(H, name='u_h(x,y)')
v = TestFunction(H)

x, y = SpatialCoordinate(mesh)
uexact = Function(H).interpolate(x * (1-x) * sin(pi * y))
uexact.rename('u(x,y) exact solution')
f = Function(H).interpolate((2 + pi**2 * x * (1-x)) * sin(pi * y))
f.rename('f(x,y)')

F = ( dot(grad(u),grad(v)) - f*v ) * dx

bdry_ids = (1, 2, 3, 4)   # four sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)

solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

# TODO:  compute difference and look at that
#udiff = Function(H).interpolate(u - uexact)
# TODO:  measure error norm
#print(errornorm(uexact, u, norm_type="L2"))

File("result.pvd").write(f,u,uexact)
