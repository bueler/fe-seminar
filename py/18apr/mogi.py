import os

os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
from firedrake import *


# Load mesh
mesh = Mesh("mogi.msh")

V = VectorFunctionSpace(mesh, "CG", 1)

# No displacement on bottom and sides (physical surface 2)
bcs = [DirichletBC(V, Constant([0, 0, 0]), 2)]

# Elastic parameters
G = 10e9
nu = 0.25
E = G * 2 * (1 + nu)

# Normal force on top of cube
P = Constant(-1e9)  # Pa

# Lame parameters from elastic constants
lam = Constant((E * nu) / ((1 + nu) * (1 - 2 * nu)))
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

# Apply pressure boundary condition inside of chamber (physical surface 1)
L = dot(P * FacetNormal(mesh), v) * ds(1)

# Solve (need to be fancier than LU decomp)
uh = Function(V, name="u")

options = {
    "ksp_type": "cg",
    "ksp_max_it": 100,
    "pc_type": "gamg",
    "pc_gamg_aggressive_square_graph": None,
    "pc_gamg_mis_k_minimum_degree_ordering": True,
    "mat_type": "aij",
    "ksp_converged_reason": None,
}

solve(a == L, uh, bcs=bcs, solver_parameters=options)

File("u.pvd").write(uh)


if True:
    # Sampling points
    xs = np.linspace(0, 5e3, 100)
    ys = np.zeros_like(xs)
    zs = np.zeros_like(xs)

    dx, dy, dz = zip(*uh.at(list(zip(xs, ys, zs))))

    # Mogi model
    a = 50   # source radius
    d = 1e3  # source depth
    dP = 1e9  # pressure change
    r = xs
    uz = (a**3)*dP*(1-nu)*d/(G*(r**2 + d**2)**(1.5))
    ur = (a**3)*dP*(1-nu)*r/(G*(r**2 + d**2)**(1.5))

    # Plot
    plt.plot(xs, dx, "k--", label="FEM Horiz.")
    plt.plot(xs, dz, "k-", label="FEM Vert.")
    plt.plot(xs, ur, "r--", label="Mogi Horiz.")
    plt.plot(xs, uz, "r-", label="Mogi Vert.")
    plt.xlabel("Radial distance (m)")
    plt.ylabel("Displacement (m)")
    plt.legend()
    plt.show()
