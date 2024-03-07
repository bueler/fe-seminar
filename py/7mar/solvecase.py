from firedrake import *
import time

defaultp = {'snes_type': 'ksponly',
            'ksp_rtol': 1.0e-10,
            'ksp_atol': 1.0e-14}
     
def solvecase(m,mesh,p):
     H = FunctionSpace(mesh,'CG',1)
     u = Function(H)
     v = TestFunction(H)

     x, y = SpatialCoordinate(mesh)
     dsqr = (x - 0.8)**2 + (y - 0.3)**2
     fsource = Function(H).interpolate(exp(-10.0*dsqr))

     F = ( dot(grad(u), grad(v)) - fsource * v ) * dx
     BCs = DirichletBC(H, Constant(0.0), (1, 2, 3, 4))

     t0 = time.time()
     solve(F == 0, u, bcs=[BCs], solver_parameters = p)
     t1 = time.time()

     dura = t1 - t0
     N = (m+1)**2
     print(f'  m = {m:5d},  N = {N:6.1e}:  {dura:7.2f} s; {1e6*dura/N:6.2f} mu s / N')
