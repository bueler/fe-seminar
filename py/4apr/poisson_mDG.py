from firedrake import *
from firedrake.output import VTKFile

# for quadrilateral mesh replace as follows:
#   mesh = UnitSquareMesh(m, m, quadrilateral=True)
#   'RT','BDM' --> 'RTCF','BDMCF'
#   'DG' --> 'DQ'

m = 20
mesh = UnitSquareMesh(m, m)

# note: some stable triangular element choices are
#   RT_k x DG_{k-1}   for k = 1,2,3,...
#   BDM_k x DG_{k-1}  for k = 1,2,3,...
# note that 'DG'-->'CG' *does* give global conservation here
k = 1
S = FunctionSpace(mesh, 'RT', k)    # or 'BDM'
H = FunctionSpace(mesh, 'DG', k-1)
W = S * H
print(f'{m} x {m} mesh mixed RT{k} x DG{k-1}:')

w = Function(W)
sigma, u = split(w)
omega, v = TestFunctions(W)

x, y = SpatialCoordinate(mesh)
f = Function(H).interpolate(5 * x + 7 * y - 1.11111)
f.rename('f')

# mixed weak form
# note Dirichlet condition on u for ids 3,4 is now "natural"
n = FacetNormal(mesh)
F = dot(sigma, omega) * dx - u * div(omega) * dx \
    + div(sigma) * v * dx - f * v * dx

# Neumann conditions on u for ids 1,2 is now Dirichlet on normal
# component of sigma = - grad(u), but we must set both components
# apparently
bc = DirichletBC(W.sub(0), as_vector([0.0,0.0]), (1,2))

solve(F == 0, w, bcs=[bc,],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu',
                           'pc_factor_mat_solver_type': 'mumps'})

sigma, u = w.subfunctions
sigma.rename('sigma')
u.rename('u')

# measure conservation success!
uint = assemble(u * dx)
fint = assemble(f * dx)
n = FacetNormal(mesh)
oflux = assemble(- dot(sigma,n) * ds)
imbalance = - oflux - fint
print(f'  u integral       = {uint:13.6e}')
print(f'  f integral       = {fint:13.6e}')
print(f'  flux out         = {oflux:13.6e}')
print(f'  imbalance        = {imbalance:13.6e}')

VTKFile("result.pvd").write(f,sigma,u)
