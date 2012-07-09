import graph3;

size(7cm);
size3(7cm);
settings.render=5; // 360 dpi

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//        2        //
//       / \       //
//      7---6      //
//     / \ / \     //
//    4---10--9    //
//   / \ / \ / \   //
//  1---5---8---3  //

triple P1 = (0, 0, 0);
triple P2 = (1, 0, 0);
triple P3 = (0, 1, 0);

triple P5 = (0, 0.32, 0.08); // (0, 0.33, 0);
triple P8 = (0, 0.68, 0.08); // (0, 0.66, 0);

triple P4 = (0.32, 0, 0.08); // (0.33, 0, 0);
triple P7 = (0.68, 0, 0.08); // (0.66, 0, 0);

triple P6 = (0.66, 0.33, 0);
triple P9 = (0.33, 0.66, 0);

triple P10 = ((P4+P5+P6+P7+P8+P9) - (P1 + P2 + P3))/3;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

path3 t = P1--P4--P7--P2--P6--P9--P3--P8--P5--P1;
path3 p1 = P4--P5--P10--P8--P9--P10--P4; 
path3 p2 = P7--P6--P10--P7;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

currentprojection = perspective(-4,2,3, center=true);

// segments
draw(P4--P1--P5, blue);
draw(P6--P2--P7, blue);
draw(P8--P3--P9, blue);

// other
draw(p1);
draw(p2);
draw(P4--P7);
draw(P5--P8);
draw(P6--P9);

dot(t,red);
dot(p1,red);
dot(p2,red);

label("$P_1$", P1, SE);
label("$P_2$", P2, E);
label("$P_3$", P3, SW);
label("$P_4$", P4, E);
label("$P_5$", P5, SE);
label("$P_6$", P6, NW);
label("$P_7$", P7, E);
label("$P_8$", P8, SE);
label("$P_9$", P9, NW);
label("$P_{10}$", P10, N);
