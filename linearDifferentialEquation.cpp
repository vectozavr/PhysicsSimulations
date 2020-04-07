//
// Created by ivan- on 08.04.2020.
//

//
// Created by ivan- on 31.03.2020.
//
#include <vector>
#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

void solve_4_equation(ComplexPlot& p1, ComplexPlot& p2, ComplexPlot& p3, ComplexPlot& p4) {

};

int main() {
    double gamma = 10;
    double E1 = 10;
    double E2 = 20;



    GnuplotPipe gp;
    gp.sendLine(R"(set yrange [-1:1])");
    gp.sendLine(R"(set xrange [0:3])");
    //gp.sendLine(gpLine);
}