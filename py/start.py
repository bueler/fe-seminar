# this getting-started code solves no differential equation problem
# it creats a mesh, computes a function f(x,y) on it,
# then it saves the result to a .pvd and plots to a figure

from firedrake import *
mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, "CG", 1)

f = Function(V)
x, y = SpatialCoordinate(mesh)
f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))
f.rename("f(x,y)")

print("writing f(x,y) to start.pvd")
File("start.pvd").write(f)

print("ploting f(x,y) in matplotlib figure")
try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

try:
  fig, axes = plt.subplots()
  from firedrake.pyplot import tripcolor
  colors = tripcolor(f, axes=axes)
  fig.colorbar(colors)
  plt.show()
except Exception as e:
  warning("Cannot plot figure. Error msg: '%s'" % e)
