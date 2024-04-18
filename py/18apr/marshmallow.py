import os

os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
from firedrake import *


# Make mesh
n = 10
mesh = UnitCubeMesh(n, n, n)

V = VectorFunctionSpace(mesh, "CG", 1)

# No displacement on bottom (physical surfaces 5)
bcs = [DirichletBC(V, Constant([0, 0, 0]), 5)]

# Elastic constants
G = 10e9  # Shear modulus (Pa)
v = 0.25  # Poisson's ratio
E = G * 2 * (1 + v)  # Young's modulus (Pa)

# Normal force on top of cube
P = Constant(-10e6)  # Pa

# Lame parameters from elastic constants
lam = Constant((E * v) / ((1 + v) * (1 - 2 * v)))
mu = Constant(G)


# Weak form
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)


def sigma(u):
    Id = Identity(mesh.geometric_dimension())
    return lam * div(u) * Id + 2 * mu * epsilon(u)


u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), grad(v)) * dx

# Normal force on top of cube (physical surface 6)
L = dot(P * FacetNormal(mesh), v) * ds(6)

# Solve
uh = Function(V, name="u")

options = {"ksp_type": "preonly", "pc_type": "lu"}

solve(a == L, uh, bcs=bcs, solver_parameters=options)

# Save deformation
File("u.pvd").write(uh)

# Stress tensor
# W = TensorFunctionSpace(mesh, "CG", 1)
# File("s.pvd").write(interpolate(sigma(uh), W))
