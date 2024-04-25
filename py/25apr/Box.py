from firedrake import *
import csv

m  = UnitSquareMesh(25,25)
mesh = ExtrudedMesh(m, 25)
x,y,z=SpatialCoordinate(mesh)
rr=(x-0.5)**2+(y-0.5)**2+(z-0.5)**2
RR=0.05

V = FunctionSpace(mesh, 'Lagrange', 1)
p = Function(V, name="p")
phi = Function(V, name="phi").interpolate(conditional(rr>RR,0,3))


u = TrialFunction(V)
v = TestFunction(V)

outfile = File("Boxout.pvd")
outfile.write(phi)

bcval = Constant(0.0)

bdry_ids = (1,2,3,4,"top","bottom")   # four sides of boundary
bc = DirichletBC(V,bcval,bdry_ids)

T = 10.
dt = 0.001
t = 0
step = 0

PHIS=[phi.dat.data[int(len(phi.dat.data)/7.)]]
TIME=[t]

while t <= T:
    step += 1
    
    phi -= dt / 2 * p
    
    solve(u * v * dx == v * p * dx + dt * inner(grad(v), grad(phi)) * dx,
        p, bcs=bc, solver_parameters={'ksp_type': 'cg',
                                    'pc_type': 'sor',
                                    'pc_sor_symmetric': True})

    phi -= dt / 2 * p

    t += dt
    
    PHIS.append(phi.dat.data[int(len(phi.dat.data)/7.)])
    TIME.append(t)
    if step % 10 == 0:
        outfile.write(phi, time=t)
      
      
  
print(PHIS)
with open('ReceiverData', 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(PHIS)

print(TIME)   
with open('ReceiverTime', 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(TIME)