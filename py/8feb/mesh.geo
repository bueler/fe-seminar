// geometry-description file for "engine block" in 2D
// to generate mesh.msh for input in block.py:
//   $ gmsh -2 mesh.geo

cl = 0.200000;  // characteristic length for triangles (elements)

// points on polygonal domain
Point(1) = {1.0,1.0,0,cl};
Point(2) = {-1.0,1.0,0,cl};
Point(3) = {-1.0,0.3,0,cl};
Point(4) = {0.0,0.3,0,cl};
Point(5) = {0.0,0.0,0,cl};
Point(6) = {-1.0,0.0,0,cl};
Point(7) = {-1.0,-1.0,0,cl};
Point(8) = {1.0,-1.0,0,cl};

// tell Gmsh about the topology and boundary of the mesh
Line(10) = {1,2};
Line(11) = {2,3};
Line(12) = {3,4};
Line(13) = {4,5};
Line(14) = {5,6};
Line(15) = {6,7};
Line(16) = {7,8};
Line(17) = {8,1};
Line Loop(20) = {10,11,12,13,14,15,16,17};
Plane Surface(30) = {20};

// indicate parts of the boundary
Physical Line(40) = {12,13,14};        // heater    (Neumann)
Physical Line(41) = {10,11,15,16,17};  // exterior  (Dirichlet)
Physical Surface(50) = {30};           // domain Omega
