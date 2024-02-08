# example showing:
#   1. reading a 2D mesh from a Gmsh file
#   2. complicated domain (more than a square ...)
#   3. Neumann boundary conditions
#   4. optional uniform mesh refinement
#   5. optional higher-order elements

from firedrake import *

# read mesh from a .msh file
#   run to get mesh.msh:  gmsh -2 mesh.geo
mesh = Mesh('mesh.msh')

# optional uniform refinement: in each refinement,
#   each triangle becomes 4 similar triangles
if False:
    refines = 2
    hierarchy = MeshHierarchy(mesh, refines)
    mesh = hierarchy[-1]

# optional higher-degree elements
degree = 1              # no higher than about 10?
H = FunctionSpace(mesh,'CG',degree)
u = Function(H, name='u(x,y)')
v = TestFunction(H)

# x, y = SpatialCoordinate(mesh)  # not used yet
# f = [put optional source term formula here]
f = Constant(0.0)

# nonhomogeneous Neumann term goes in weak form
# gN = [put optional boundary heating formula here]
gN = Constant(-10.0)
F = ( dot(grad(u),grad(v)) - f * v ) * dx + gN * v * ds

# explicit: part of the boundary is Dirichlet
exterior = (41,)  # see mesh.geo for this boundary id
# gD = [put optional boundary temperature formula here]
gD = Constant(0.0)
BCs = DirichletBC(H, gD, exterior)

# solve, and write solution, as before
solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})
File("result.pvd").write(u)
