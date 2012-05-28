import graph3;

size(7cm);
size3(7cm);
settings.render=5; // 360 dpi

currentprojection = perspective(-4,2,1, center=true);

draw((0,0,0)--(-1,0,0)--(0,-1,0)--cycle);

// theta_x
draw((0,0,0)--(-0.2,0,0), blue, EndArrow3);
label("$\theta_x$", (-0.2,0,0), SW, blue);

// theta_y
draw((0,0,0)--(0,-0.2,0), blue, EndArrow3);
label("$\theta_y$", (0,-0.2,0), NW, blue);

// u_z
draw((0,0,0)--(0,0,0.2), green, EndArrow3);
label("$u_z$", (0,0,0.2), W, green);
