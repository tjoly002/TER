// Gmsh project created on Tue Dec  5 11:54:07 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {3.5, 0, 0, 1.0};
//+
Point(3) = {3.5, 40, 0, 1.0};
//+
Point(4) = {0, 40, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("NordOuestEst", 10) = {4, 3, 2};
//+
Physical Curve("Sud", 11) = {1};
//+
Physical Surface("Tinit", 100) = {1};
