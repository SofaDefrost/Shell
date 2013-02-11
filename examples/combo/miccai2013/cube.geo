lc1 = 0.018;
lc2 = 0.018;

// rectangle

pxl=0;
pxh=1;
pyl=0;
pyh=1;
pzl=0;
pzh=1;

Point(1) = {pxl,pyl,pzl,lc1};
Point(2) = {pxh,pyl,pzl,lc1};
Point(3) = {pxh,pyh,pzl,lc1};
Point(4) = {pxl,pyh,pzl,lc1};

Point(11) = {pxl,pyl,pzh,lc2};
Point(12) = {pxh,pyl,pzh,lc2};
Point(13) = {pxh,pyh,pzh,lc2};
Point(14) = {pxl,pyh,pzh,lc2};

// bottom face
Line(101) = {4,3};
Line(102) = {3,2};
Line(103) = {2,1};
Line(104) = {1,4};
Line Loop(200) = {102,103,104,101};
Plane Surface(1200) = {200};

// top face
Line(111) = {14,13};
Line(112) = {13,12};
Line(113) = {12,11};
Line(114) = {11,14};
Line Loop(201) = {112,113,114,111};
Plane Surface(1201) = {201};

// vertical
Line(131) = {11,1};
Line(132) = {12,2};
Line(133) = {13,3};
Line(134) = {14,4};
Line Loop(300) = {103, -131, -113, 132};
Plane Surface(1300) = {300};

Line Loop(301) = {104, -134, -114, 131};
Plane Surface(1301) = {301};

Line Loop(302) = {101, -133, -111, 134};
Plane Surface(1302) = {302};

Line Loop(303) = {102, -132, -112, 133};
Plane Surface(1303) = {303};

// bottom cube:
Surface Loop(10000) = {1301,1200,1300,1303,1302,1201};
Volume(20000) = {10000};
Pysical Volume(30000) = {20000};
