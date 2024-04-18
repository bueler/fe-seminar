import gmsh
import rasterio as rio
import numpy as np
from tqdm import tqdm
import sys


# Extrude mesh up from plane to match topography?


# GMSH boilerplate
gmsh.initialize()
gmsh.clear()
lc = 50
dep = 5000

# Conduit

cx, cy = (32270.0, 1040901.8)  # x,y
cr = 25  # radius
ct = 1000  # top
cb = -1000  # bottom

# Model name
gmsh.model.add("topo")

fd = rio.open("./augustine_ifsar_clip_100m.tif")
dem = fd.read(1)
n, m = dem.shape
gt = fd.transform
lft, bot, rgt, top = fd.bounds
fd.close()

# Make grid of points
points = np.empty(dem.shape, dtype=int)

for i in tqdm(range(0, dem.shape[0])):
    for j in range(0, dem.shape[1]):
        x, y = gt*(j, i)
        z = dem[i, j]
        points[i, j] = gmsh.model.geo.addPoint(x, y, z, lc)

# Make horiz, vert, diag lines
horiz = np.empty((points.shape[0], points.shape[1]-1), dtype=int)
vert = np.empty((points.shape[0]-1, points.shape[1]), dtype=int)
diag = np.empty((points.shape[0]-1, points.shape[1]-1), dtype=int)

for i in tqdm(range(0, points.shape[0])):
    for j in range(0, points.shape[1]-1):
        horiz[i,j] = gmsh.model.geo.addLine(points[i,j], points[i,j+1])

for i in tqdm(range(0, points.shape[0]-1)):
    for j in range(0, points.shape[1]):
        vert[i,j] = gmsh.model.geo.addLine(points[i,j], points[i+1,j])

for i in tqdm(range(0, points.shape[0]-1)):
    for j in range(0, points.shape[1]-1):
        diag[i,j] = gmsh.model.geo.addLine(points[i+1,j+1], points[i,j])

srfs = []

# Connect into surfaces
for i in tqdm(range(0, points.shape[0]-1)):
    for j in range(0, points.shape[1]-1):
        cl0 = gmsh.model.geo.addCurveLoop([horiz[i,j], vert[i,j+1], diag[i,j]])
        cl1 = gmsh.model.geo.addCurveLoop([vert[i,j], horiz[i+1,j], diag[i,j]])
        srfs.append(gmsh.model.geo.addPlaneSurface([cl0]))
        srfs.append(gmsh.model.geo.addPlaneSurface([cl1]))

# Add sides and bottom
x, y = gt*(0, 0)
p0 = gmsh.model.geo.addPoint(x, y, -dep, lc)

x, y = gt*(m-1, 0)
p1 = gmsh.model.geo.addPoint(x, y, -dep, lc)

x, y = gt*(m-1, n-1)
p2 = gmsh.model.geo.addPoint(x, y, -dep, lc)

x, y = gt*(0, n-1)
p3 = gmsh.model.geo.addPoint(x, y, -dep, lc)

l0 = gmsh.model.geo.addLine(p0, p1)
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p0)

l4 = gmsh.model.geo.addLine(p0, points[0, 0])
l5 = gmsh.model.geo.addLine(p1, points[0, m-1])
l6 = gmsh.model.geo.addLine(p2, points[n-1, m-1])
l7 = gmsh.model.geo.addLine(p3, points[n-1, 0])

# Floor
cl0 = gmsh.model.geo.addCurveLoop([l0, l1, l2, l3])
ps0 = gmsh.model.geo.addPlaneSurface([cl0])

# Left
cl1 = gmsh.model.geo.addCurveLoop([l4] + list(vert[:, 0]) + [-l7, l3])
ps1 = gmsh.model.geo.addPlaneSurface([cl1])

# Right
cl2 = gmsh.model.geo.addCurveLoop([l5] + list(vert[:, m-1]) + [-l6, -l1])
ps2 = gmsh.model.geo.addPlaneSurface([cl2])

# Bottom
cl3 = gmsh.model.geo.addCurveLoop([l7] + list(horiz[n-1, :]) + [-l6, l2])
ps3 = gmsh.model.geo.addPlaneSurface([cl3])

# Top
cl4 = gmsh.model.geo.addCurveLoop([l4] + list(horiz[0, :]) + [-l5, -l0])
ps4 = gmsh.model.geo.addPlaneSurface([cl4])

# Make volume
sl0 = gmsh.model.geo.addSurfaceLoop([ps0, ps1, ps2, ps3, ps4] + srfs)
#v0 = gmsh.model.geo.addVolume([sl0])

# Make conduit top
p4 = gmsh.model.geo.addPoint(cx, cy, ct, lc)
p5 = gmsh.model.geo.addPoint(cx-cr, cy, ct, lc)
p6 = gmsh.model.geo.addPoint(cx, cy+cr, ct, lc)
p7 = gmsh.model.geo.addPoint(cx+cr, cy, ct, lc)
p8 = gmsh.model.geo.addPoint(cx, cy-cr, ct, lc)

l8 = gmsh.model.geo.addCircleArc(p5, p4, p6)
l9 = gmsh.model.geo.addCircleArc(p6, p4, p7)
l10 = gmsh.model.geo.addCircleArc(p7, p4, p8)
l11 = gmsh.model.geo.addCircleArc(p8, p4, p5)

cl5 = gmsh.model.geo.addCurveLoop([l8, l9, l10, l11])
ps5 = gmsh.model.geo.addPlaneSurface([cl5])

# Extrude sides
ex0 = gmsh.model.geo.extrude([(1, l8)], 0, 0, cb-ct)
ex1 = gmsh.model.geo.extrude([(1, l9)], 0, 0, cb-ct)
ex2 = gmsh.model.geo.extrude([(1, l10)], 0, 0, cb-ct)
ex3 = gmsh.model.geo.extrude([(1, l11)], 0, 0, cb-ct)

# Bottom
cl6 = gmsh.model.geo.addCurveLoop([ex0[0][1], ex1[0][1], ex2[0][1], ex3[0][1]])
ps6 = gmsh.model.geo.addPlaneSurface([cl6])

# Conduit surface loop
sl1 = gmsh.model.geo.addSurfaceLoop([ps5, ex0[1][1], ps6, ex1[1][1], ex2[1][1], ex3[1][1]])

v0 = gmsh.model.geo.addVolume([sl0, sl1])

gmsh.model.geo.synchronize()

# Set mesh size as a function of distance from pressure source
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.model.mesh.field.add("MathEval", 1)
gmsh.model.mesh.field.setString(1, "F","0.1*sqrt((x-%f)^2 + (y-%f)^2) + 10" % (cx, cy))

gmsh.model.mesh.field.setAsBackgroundMesh(1)


# Add physical groups
gmsh.model.addPhysicalGroup(2, [ps0, ps1, ps2, ps3, ps4], tag=1)  # Dirichlet condition

# With top and bottom
#gmsh.model.addPhysicalGroup(2, [ps5, ex0[1][1], ps6, ex1[1][1], ex2[1][1], ex3[1][1]], tag=2) # Pressure

# Without top and bottom
gmsh.model.addPhysicalGroup(2, [ex0[1][1], ex1[1][1], ex2[1][1], ex3[1][1]], tag=2) # Pressure

gmsh.model.addPhysicalGroup(3, [v0], tag=3)  # volume

gmsh.model.mesh.generate(3)

gmsh.write("conduit.msh")

gmsh.finalize()