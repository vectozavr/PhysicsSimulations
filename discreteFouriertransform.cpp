//
// Created by ivan- on 11.02.2020.
//
#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;

int main() {
    ComplexPlot data;
    ComplexPlot data2;
    ComplexPlot transform;

    int quantity = 1000;
    for(int i = -quantity; i < quantity; i++) {
        if(i < -quantity / 10 || i > quantity / 10)
            data.push((double)i/quantity*10, {0, 0});
        else
            data.push((double)i/quantity*10, {1, 0});
    }
    //for(int i = -quantity; i < quantity; i++) {
    //    double x = (double)i/quantity * 2 * PI;
    //    data.push_back({x, sin(x) + sin(10*x)/2});
    //}

    //addNoise(data);
    fourierTransform(data, transform);
    inverseFourierTransform(transform, data2);

    saveVectorPoint2DToFile(data.abs(), "data.dat");
    saveVectorPoint2DToFile(data2.abs(), "transform.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");
    gp.sendLine(R"(set multiplot layout 1,2)");
    gp.sendLine(R"(plot "data.dat" with lines)");
    gp.sendLine(R"(plot "transform.dat" with lines)");

    return 0;

}