lc1 = 0.1;
lc2 = 0.015;
lc3 = 0.015;

// rectangle

pxl=0;
pxh=1;
pyl=0;
pyh=1;
pzl=0;
pzh=0.92;
pzhh=1.00;

Point(1) = {pxl,pyl,pzl,lc1};
Point(2) = {pxh,pyl,pzl,lc1};
Point(3) = {pxh,pyh,pzl,lc1};
Point(4) = {pxl,pyh,pzl,lc1};

Point(11) = {pxl,pyl,pzh,lc2};
Point(12) = {pxh,pyl,pzh,lc2};
Point(13) = {pxh,pyh,pzh,lc2};
Point(14) = {pxl,pyh,pzh,lc2};

Point(21) = {pxl,pyl,pzhh,lc3};
Point(22) = {pxh,pyl,pzhh,lc3};
Point(23) = {pxh,pyh,pzhh,lc3};
Point(24) = {pxl,pyh,pzhh,lc3};

// bottom face
Line(101) = {4,3};
Line(102) = {3,2};
Line(103) = {2,1};
Line(104) = {1,4};
Line Loop(200) = {102,103,104,101};
Plane Surface(1200) = {200};

//middle face
Line(111) = {14,13};
Line(112) = {13,12};
Line(113) = {12,11};
Line(114) = {11,14};
Line Loop(201) = {112,113,114,111};
Plane Surface(1201) = {201};

//top face
Line(121) = {24,23};
Line(122) = {23,22};
Line(123) = {22,21};
Line(124) = {21,24};
Line Loop(202) = {122,123,124,121};
Plane Surface(1202) = {202};

//bottom cube vertical
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

// top cube vertical
Line(141) = {21,11};
Line(142) = {22,12};
Line(143) = {23,13};
Line(144) = {24,14};

Line Loop(400) = {113, -141, -123, 142};
Plane Surface(1400) = {400};

Line Loop(401) = {114, -144, -124, 141};
Plane Surface(1401) = {401};

Line Loop(402) = {111, -143, -121, 144};
Plane Surface(1402) = {402};

Line Loop(403) = {112, -142, -122, 143};
Plane Surface(1403) = {403};


// bottom cube:
Surface Loop(10000) = {1301,1200,1300,1303,1302,1201};
Volume(20000) = {10000};
Pysical Volume(30000) = {20000};

Surface Loop(10001) = {1401,1201,1400,1403,1402,1202};
Volume(20001) = {10001};
Pysical Volume(30001) = {20001};
