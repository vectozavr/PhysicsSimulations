//
// Created by Иван Ильин on 10.08.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

int main() {
    ComplexPlot p_R;

    GnuplotPipe gp;

    //gp.sendLine(R"(set xrange [-15:15])");
    //gp.sendLine(R"(set yrange [-15:15])");
    gp.sendLine(R"(plot "log.txt" with lines)");
}