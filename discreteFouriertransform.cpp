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
    //for(int i = -quantity; i < quantity; i++) {
    //    if(i < -quantity / 10 || i > quantity / 10)
    //        data.push((double)i/quantity*10, {0, 0});
    //    else
    //        data.push((double)i/quantity*10, {1, 0});
    //}
    for(int i = -quantity; i < quantity; i++) {
        double x = (double)i/quantity * 10;
        if(i < 0)
            data.push(x, {sin(x), 0});
        else
            data.push(x, {sin(-x), 0});
    }

    //addNoise(data);
    fourierTransform(data, transform);
    //inverseFourierTransform(transform, data2);

    saveVectorPoint2DToFile(data.real(), "data.dat");
    saveVectorPoint2DToFile(transform.abs(), "transform.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");
    //gp.sendLine(R"(set multiplot layout 1,2)");
    //gp.sendLine(R"(plot "data.dat" with lines)");
    gp.sendLine(R"(plot "transform.dat" with lines)");

    return 0;

}