// Gmsh project created on Wed Nov 29 12:51:52 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 35, 0, 1.0};
//+
Point(3) = {0.5, 35, 0, 1.0};
//+
Point(4) = {0.5, 0, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Sud", 11) = {2};
//+
Physical Curve("NordOuestEst", 10) += {3, 1, 4};
//+
Physical Surface("Tinit", 100) = {1};
//+
Transfinite Surface {1} Right;
//+
Recombine Surface {1};
//+
Transfinite Surface {1};
