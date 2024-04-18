import os

os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
from firedrake import *


# Load mesh
mesh = Mesh("conduit.msh")

V = VectorFunctionSpace(mesh, "CG", 1)

# No displacement on bottom and sides (physical surface 1)
bcs = [DirichletBC(V, Constant([0, 0, 0]), 1)]

# Elastic parameters
G = 10e9
v = 0.25
E = G*2*(1+v)

# Pressure in chamber
P = Constant(-10e7)

# Lame parameters
lam = Constant((E*v)/((1+v)*(1-2*v)))
mu = Constant(G)
Id = Identity(mesh.geometric_dimension())


def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)


def sigma(u):
    return lam * div(u) * Id + 2 * mu * epsilon(u)


u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), grad(v)) * dx

# Apply pressure boundary condition inside of chamber (physical surface 2)
L = dot(P*FacetNormal(mesh), v) * ds(2)

# Solve
uh = Function(V, name="u")

options = {"ksp_type": "cg", 
           "ksp_max_it": 100, 
           "pc_type": "gamg",
           "pc_gamg_aggressive_square_graph": None,
           "pc_gamg_mis_k_minimum_degree_ordering": True,
           "mat_type": "aij",
           "ksp_converged_reason": None}

solve(a == L, uh, bcs=bcs, solver_parameters=options)

File("u.pvd").write(uh)

File("tilt.pvd").write(interpolate(grad(uh.sub(2)), V))
