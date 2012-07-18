import graph3;

size(7cm);
size3(7cm);
settings.render=5; // 360 dpi

triple P1=(0,0,0); // don't change this one
triple P2=(0.5,0.866,0);
triple P3=(1,0,0);

currentprojection = perspective((1,-4,4), center=true, target=(P1+P2+P3)/3);

draw(P1--P2--P3--cycle);

// theta_x
draw(P1--(0.3,0,0), blue, EndArrow3);
label("$\theta_x$", (0.3,0,0), SW, blue);

// theta_y
draw(P1--(0,0.3,0), blue, EndArrow3);
label("$\theta_y$", (0,0.3,0), NW, blue);

// u_z
draw(P1--(0,0,0.3), green, EndArrow3);
label("$u_z$", (0,0,0.3), W, green);
