# This example solves a linear system of PDEs, i.e. a vector PDE:
#   - v_xx - v_yy               = f(x,y)
#   - v           - u_xx - u_yy = 0
# for  w = [v(x,y); u(x,y)].  The discretized system produces a block-
# structured matrix problem
#   [  L   0 ] = [ f ]
#   [ -I   L ] = [ 0 ]
# where L is the usual discretization of the Laplacian.  This system
# models an elastic plate, and u(x,y) satisfies the scalar fourth-order
# "pure bending" plate equation:
#   u_xxxx + 2 u_xxyy + u_yyyy = 0
# (See: https://en.wikipedia.org/wiki/Plate_theory)
# This code demonstrates basic "fieldsplit" preconditioning options to
# solve the block-structured system by an additive multigrid method.
# Documentation:  slides/21mar-plate.pdf

from firedrake import *
from firedrake.output import VTKFile  # in older Firedrake versions: delete this line

M = 40    # double M several times to see solver scaling
          # at M=640 the multigrid solver is distinctly faster

# mesh and 2 component function space
mesh = UnitSquareMesh(M,M)
H = FunctionSpace(mesh, 'CG', 1)
W = H * H            # w = [v; u] is in this space
w = Function(W, name='w_h(x,y) vector solution')
v, u = split(w)
r, s = TestFunctions(W)

# load
# x, y = SpatialCoordinate(mesh)
f = Constant(1.0)  # FIXME set up exact solution from Chapter 7 example

# weak form:  F == 0  for all r,s
F = ( dot(grad(v), grad(r)) - f * r \
      - v * s + dot(grad(u), grad(s)) ) * dx

bdry_ids = (1, 2, 3, 4)   # four sides of boundary
BCs = DirichletBC(H, Constant(0.0), bdry_ids)
BCs = [ DirichletBC(W.sub(0), Constant(0.0), bdry_ids),
        DirichletBC(W.sub(1), Constant(0.0), bdry_ids) ]

# any solver can use these options
par = {'snes_type': 'ksponly',
    'ksp_type': 'gmres',
    #'ksp_view': None,
    #'ksp_monitor': None,
    'ksp_converged_reason': None}

# choose between solvers by par.update(ONE OF THESE):
#   monolithic direct
directLU = {'pc_type': 'lu'}
#   fieldsplit direct
FSmulLU = {'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'multiplicative', # 'additive' will need 2 iterations
    'fieldsplit_0_ksp_type': 'preonly',
    'fieldsplit_0_pc_type': 'lu',
    'fieldsplit_1_ksp_type': 'preonly',
    'fieldsplit_1_pc_type': 'lu'}
#   fieldsplit algebraic multigrid  (scalable multigrid-preconditioned iterative)
FSaddMG = {'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'additive', # 'multiplicative' will need fewer iterations
    'fieldsplit_0_ksp_type': 'preonly',
    'fieldsplit_0_pc_type': 'gamg',
    'fieldsplit_1_ksp_type': 'preonly',
    'fieldsplit_1_pc_type': 'gamg'}
par.update(directLU)
#par.update(FSmulLU)
#par.update(FSaddMG)

solve(F == 0, w, bcs=BCs,
      options_prefix='s', # add -s_xxx_yyy options at runtime
      solver_parameters = par)

v = w.subfunctions[0]
v.rename('v_h(x,y) moment')
u = w.subfunctions[1]
u.rename('u_h(x,y) displacement')
normv = norm(v, norm_type='L2')
normu = norm(u, norm_type='L2')
print(f'{M:d} x {M:d} mesh: |v| = {normv:.6f}, |u| = {normu:.6f}')

print('writing solution v, u to result.pvd ...')
VTKFile("result.pvd").write(v, u)  # older versions of Firedrake use "File()"
