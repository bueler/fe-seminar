from firedrake import *
from solvecase import solvecase, defaultp

levs = [4,5,6,7,8,9]  # coarse mesh is 4 x 4
#levs = [4,5,6,7,8,9,10]  # WARNING memory?
p = defaultp.copy()
p.update({'ksp_type': 'cg',
          'pc_type': 'mg',
          'mg_coarse_ksp_type': 'preonly',
          'mg_coarse_pc_type': 'lu',
          'mg_levels_ksp_type': 'richardson',
          'mg_levels_ksp_max_it': 2,
          'mg_levels_pc_type': 'icc',
          #'ksp_monitor': None,
          #'ksp_view': None,
          'ksp_converged_reason': None})

print('solve time for m x m meshes with N dofs:')
for k in levs:
    coarse = UnitSquareMesh(4,4)
    hierarchy = MeshHierarchy(coarse, k)
    mesh = hierarchy[-1]
    m = 4 * 2**k
    solvecase(m, mesh, p)
