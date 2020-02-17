//
// Created by ivan- on 17.02.2020.
//

#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;

int main() {
    std::vector<Point2D> data;
    std::vector<Point2D> data2;
    std::vector<Point2D> transform;

    int quantity = 1000;
    for(int i = -quantity; i < quantity; i++) {
        if(i < -quantity / 10 || i > quantity / 10)
            data.push_back({(double)i/quantity*10, 0});
        else
            data.push_back({(double)i/quantity*10, 1});
    }
    //for(int i = -quantity; i < quantity; i++) {
    //    double x = (double)i/quantity * 2 * PI;
    //    data.push_back({x, sin(x)});
    //}

    //addNoise(data);
    //fourierTransform(data, transform);
    //inverseFourierTransform(transform, data2);

    saveVectorPoint2DToFile(data, "data.dat");
    saveVectorPoint2DToFile(data2, "transform.dat");

    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 1,2)");
    gp.sendLine(R"(plot "data.dat" with lines)");
    gp.sendLine(R"(plot "transform.dat" with lines)");

    return 0;
}