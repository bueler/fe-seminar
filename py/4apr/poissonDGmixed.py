from firedrake import *
from firedrake.output import VTKFile

m = 20
mesh = UnitSquareMesh(m, m)

k = 1
print(f'{m} x {m} mesh mixed RT{k} x DG{k-1}:')
# note: some stable choices are
#   RT_k x DG_{k-1}   for k = 1,2,3,...
#   BDM_k x DG_{k-1}  for k = 1,2,3,...
S = FunctionSpace(mesh, 'RT', k)
H = FunctionSpace(mesh, 'DG', k-1)
W = S * H

w = Function(W)
sigma, u = split(w)
omega, v = TestFunctions(W)

x, y = SpatialCoordinate(mesh)
f = Function(H).interpolate(5 * x + 7 * y - 1.11111)
f.rename('f')

# mixed weak form
# Dirichlet condition on u for ids 3,4 is now "natural"
n = FacetNormal(mesh)
F = dot(sigma,omega) * dx - u * div(omega) * dx + u * dot(omega,n) * ds((1,2)) \
    + div(sigma) * v * dx - f * v * dx

# Neumann conditions on u for ids 1,2 is now Dirichlet on sigma = - grad(u)
bc = DirichletBC(W.sub(0), as_vector([0.0,0.0]), (1,2))

solve(F == 0, w, bcs=[bc,],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu',
                           'pc_factor_mat_solver_type': 'mumps'})

sigma, u = w.subfunctions
sigma.rename('sigma')
u.rename('u')

# measure conservation failure (is success!)
uint = assemble(u * dx)
fint = assemble(f * dx)
n = FacetNormal(mesh)
oflux = assemble(- dot(sigma,n) * ds)
imbalance = - oflux - fint
print(f'  u integral       = {uint:10.3e}')
print(f'  f integral       = {fint:10.3e}')
print(f'  flux out         = {oflux:10.3e}')
print(f'  imbalance        = {imbalance:10.3e}')

VTKFile("result.pvd").write(f,sigma,u)
