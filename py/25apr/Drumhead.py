from firedrake import *

#Simple Square Mesh
mesh  = UnitSquareMesh(25,25)

#generating coordinates for step function pulse
x,y=SpatialCoordinate(mesh)

#conditon will be a circlular pulse of radius 0.05
rr=(x-0.5)**2+(y-0.5)**2
RR=0.05

V = FunctionSpace(mesh, 'Lagrange', 1)

#need both p and phi for verlet time stepping method
p = Function(V, name="p") 
phi = Function(V, name="phi").interpolate(conditional(rr>RR,0,3)) #initial condition 

u = TrialFunction(V)
v = TestFunction(V)

outfile = File("Drumheadout.pvd")
outfile.write(phi)

#0 diriclet BC
bcval = Constant(0.0)

bdry_ids = (1,2,3,4)   # four sides of boundary
bc = DirichletBC(V,bcval, bdry_ids)

#setting time loop parameters
T = 10.
dt = 0.001
t = 0
step = 0

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
    if step % 10 == 0:
        outfile.write(phi, time=t) #don't need to save every timestep to animate