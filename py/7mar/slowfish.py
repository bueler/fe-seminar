from firedrake import *
from solvecase import solvecase, defaultp

# uncomment lines below to save the internal matrix A to A.m

M = [64,128,256,512,1024]
#M = [8,]
p = defaultp.copy()
p.update({#'ksp_view_mat': ':A.m:ascii_matlab',
          #'ksp_view': None,
          'ksp_type': 'preonly',
          'pc_type': 'lu'})

print('solve time for m x m meshes with N dofs:')
for m in M:
    mesh = UnitSquareMesh(m,m)
    solvecase(m,mesh,p)
