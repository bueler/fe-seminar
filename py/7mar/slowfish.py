from firedrake import *
from solvecase import solvecase, defaultp

M = [64,128,256,512,1024]
p = defaultp.copy()
p.update({'ksp_type': 'preonly',
          'pc_type': 'lu'})

print('solve time for m x m meshes with N dofs:')
for m in M:
    mesh = UnitSquareMesh(m,m)
    solvecase(m,mesh,p)
