//
// Created by ivan- on 17.02.2020.
//

#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;

int main() {
    ComplexPlot data1;  // step
    ComplexPlot data2;  // step + noise
    ComplexPlot data3;  // sin10
    ComplexPlot data4;  // sin20 + sin10
    ComplexPlot data5;  // sin10 + noise
    ComplexPlot data6;  // sin9
    ComplexPlot data7;  // shift phase

    ComplexPlot conv1;
    ComplexPlot conv2;
    ComplexPlot conv3;
    ComplexPlot conv4;
    ComplexPlot conv5;
    ComplexPlot conv6;
    ComplexPlot conv7;
    ComplexPlot transformCheck;


    int quantity = 500;
    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;

        if(i < quantity / 2 || i > 6 * quantity / 10)
            data1.push(x, {0, 0});
        else
            data1.push(x, {1, 0});

        data3.push(x, {sin(10*x), 0});
        data4.push(x, {sin(20*x) + sin(10*x), 0});

        data6.push(x, {sin(9*x), 0});

        if(i < quantity / 2)
            data7.push(x, {sin(5*x), 0});
        else
            data7.push(x, {sin(-5*x), 0});
    }
    data2 = data1;
    data5 = data3;
    addNoise(data2);
    addNoise(data5);

    convolution(data1, data1, conv1);
    convolution(data2, data2, conv2);
    convolution(data1, data3, conv3);
    convolution(data3, data4, conv4);

    saveVectorPoint2DToFile(data1.real(), "data1.dat");
    saveVectorPoint2DToFile(data2.real(), "data2.dat");
    saveVectorPoint2DToFile(data3.real(), "data3.dat");
    saveVectorPoint2DToFile(data4.real(), "data4.dat");

    saveVectorPoint2DToFile(conv1.real(), "conv1.dat", conv1.size()/2);
    saveVectorPoint2DToFile(conv2.real(), "conv2.dat", conv2.size()/2);
    saveVectorPoint2DToFile(conv3.real(), "conv3.dat", conv3.size()/2);
    saveVectorPoint2DToFile(conv4.real(), "conv4.dat", conv4.size()/2);

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");

    gp.sendLine(R"(set multiplot layout 4,3)");
    gp.sendLine(R"(unset key)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "conv1.dat" with lines)");

    gp.sendLine(R"(plot "data2.dat" with lines)");
    gp.sendLine(R"(plot "data2.dat" with lines)");
    gp.sendLine(R"(plot "conv2.dat" with lines)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "conv3.dat" with lines)");

    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "data4.dat" with lines)");
    gp.sendLine(R"(plot "conv4.dat" with lines)");

    return 0;
}