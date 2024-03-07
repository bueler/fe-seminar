from firedrake import *
from solvecase import solvecase, defaultp

M = [100,200,400,800,1200]
p = defaultp.copy()
p.update({'ksp_type': 'cg',
          'pc_type': 'icc',
          #'ksp_monitor': None,
          'ksp_converged_reason': None,
          'ksp_rtol': 1.0e-12})

print('solve time for m x m mesh')
for m in M:
     mesh = UnitSquareMesh(m,m)
     solvecase(m,mesh,p)
