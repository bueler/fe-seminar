from firedrake import *
import csv #for writing timeseries output

#making an extruded, 3D mesh
m = IcosahedralSphereMesh(radius=10, refinement_level=5)
mesh = ExtrudedMesh(m, 25, layer_height=0.1, extrusion_type='radial')


#conditon will be a circlular pulse of radius 0.02
#located near the bottom of the domain
x,y,z=SpatialCoordinate(mesh)
rr=(x)**2+(y)**2+(z-10.25)**2
RR=0.02

#need both p and phi for verlet time stepping method
V = FunctionSpace(mesh, 'Lagrange', 1)
p = Function(V, name="p")
phi = Function(V, name="phi").interpolate(conditional(rr>RR,0,3))

u = TrialFunction(V)
v = TestFunction(V)

outfile = File("EIWaveout.pvd")
outfile.write(phi)

#0 diriclet BC
bcval = Constant(0.0)

# For this mesh, there is only the top and bottom
bdry_ids = ("top","bottom")
bc = DirichletBC(V,bcval,bdry_ids)

#setting time loop parameters
T = 30
dt = 0.001
t = 0
step = 0

#Saving timeseries
PHIS=[phi.dat.data[int(len(phi.dat.data)/7.)]]
TIME=[t]

#time loop
while t <= T:
    step += 1
    
    #step one of verlet method
    phi -= dt / 2 * p
    
    #step two of verlet method
    #invert mass matrix using linear solver
    solve(u * v * dx == v * p * dx + dt * inner(grad(v), grad(phi)) * dx,
        p, bcs=bc, solver_parameters={'ksp_type': 'cg',
                                    'pc_type': 'sor',
                                    'pc_sor_symmetric': True})

    #step three of verlet method
    phi -= dt / 2 * p
    t += dt
    
    #Saving timeseries
    PHIS.append(phi.dat.data[int(len(phi.dat.data)/7.)])
    TIME.append(t)
    if step % 10 == 0:
        outfile.write(phi, time=t)

#Writing timeseries
print(PHIS)
with open('ReceiverData', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(PHIS)
 
print(TIME)   
with open('ReceiverTime', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(TIME)