# The commented lines at bottom are from the last slide for
# the Thursday 18 January seminar.  I will edit them as
# to solve a particular Poisson problem.  See poisson.py
# for the cleaned-up and final version of the demo.

# Recall that the Poisson equation
#   - div grad u = f
# with zero boundary conditions becomes the weak form
#   int_Omega grad u . grad v - f v dx = 0
# for all v in H_0^1(Omega).

from firedrake import *

#mesh = UnitSquareMesh(10,10)
#H = FunctionSpace(mesh,'CG',1)
#u = Function(H)
#v = TestFunction(H)
#f = Function(H).interpolate(...)
#F = ( dot(grad(u),grad(v)) - f*v ) * dx
#BCs = ...
#solve(F == 0, u, bcs = [BCs], ...)
