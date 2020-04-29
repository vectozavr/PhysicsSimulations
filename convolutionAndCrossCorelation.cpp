//
// Created by ivan- on 17.02.2020.
//

#include "gnuplot.h"
#include "vemath.h"
#include "fftw/fftw3.h"

using namespace vemath;

int main() {
    ComplexPlot data1;  // step
    ComplexPlot data2;  // step + noise
    ComplexPlot data3;  // sin10
    ComplexPlot data4;  // sin20 + sin10
    ComplexPlot data5;  // sin10 + noise
    ComplexPlot data6;  // sin9
    ComplexPlot data7;  // shift phase
    ComplexPlot data9;  // step2
    ComplexPlot data10; // sin20
    ComplexPlot data11; // noise
    ComplexPlot data12; // sin
    ComplexPlot data13; // sin3

    ComplexPlot conv1;
    ComplexPlot conv2;
    ComplexPlot conv3;
    ComplexPlot conv4;
    ComplexPlot conv5;
    ComplexPlot conv6;
    ComplexPlot conv7;
    ComplexPlot conv8;
    ComplexPlot conv9;
    ComplexPlot conv10;
    ComplexPlot transformCheck;


    int quantity = 500;
    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;

        data11.push(x, {0, 0});

        if(i < quantity / 2 || i > 6 * quantity / 10)
            data1.push(x, {0, 0});
        else
            data1.push(x, {1, 0});

        if(i < quantity / 2 || i > 8 * quantity / 10)
            data9.push(x, {0, 0});
        else
            data9.push(x, {1, 0});

        data3.push(x, {sin(10*x), 0});
        data4.push(x, {sin(20*x) + sin(10*x), 0});

        data6.push(x, {sin(9*x), 0});
        data10.push(x, {sin(20*x), 0});

        data12.push(x, {sin(x), 0});
        data13.push(x, {sin(3*x), 0});

        if(i < quantity / 2)
            data7.push(x, {sin(5*x), 0});
        else
            data7.push(x, {sin(-5*x), 0});
    }
    data2 = data1;
    data5 = data3;
    addNoise(data2);
    addNoise(data5);
    addNoise(data11);

    convolution(data1, data1, conv1);
    convolution(data2, data2, conv2); // step + noise
    //crossCorrelation(data2, data2, conv2); // step + step cross correlation
    convolution(data1, data3, conv3);
    convolution(data3, data3, conv4);
    convolution(data1, data9, conv5);
    convolution(data3, data10, conv6);
    convolution(data11, data11, conv7);
    convolution(data1, data2, conv8);
    convolution(data5, data1, conv9);
    convolution(data12, data13, conv10);

    saveVectorPoint2DToFile(data1.real(), "data1.dat");
    saveVectorPoint2DToFile(data2.real(), "data2.dat");
    saveVectorPoint2DToFile(data3.real(), "data3.dat");
    saveVectorPoint2DToFile(data4.real(), "data4.dat");
    saveVectorPoint2DToFile(data9.real(), "data9.dat");
    saveVectorPoint2DToFile(data10.real(), "data10.dat");
    saveVectorPoint2DToFile(data11.real(), "data11.dat");
    saveVectorPoint2DToFile(data12.real(), "data12.dat");
    saveVectorPoint2DToFile(data13.real(), "data13.dat");

    saveVectorPoint2DToFile(conv1.real(), "conv1.dat");
    saveVectorPoint2DToFile(conv2.real(), "conv2.dat");
    saveVectorPoint2DToFile(conv3.real(), "conv3.dat");
    saveVectorPoint2DToFile(conv4.real(), "conv4.dat");
    saveVectorPoint2DToFile(conv5.real(), "conv5.dat");
    saveVectorPoint2DToFile(conv6.real(), "conv6.dat");
    saveVectorPoint2DToFile(conv7.real(), "conv7.dat");
    saveVectorPoint2DToFile(conv8.real(), "conv8.dat");
    saveVectorPoint2DToFile(conv9.real(), "conv9.dat");
    saveVectorPoint2DToFile(conv10.real(), "conv10.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");

    gp.sendLine(R"(set multiplot layout 5,4)");
    gp.sendLine(R"(unset key)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "conv1.dat" with lines)");
    gp.sendLine(R"(plot "conv7.dat" with lines)");

    gp.sendLine(R"(plot "data2.dat" with lines)");
    gp.sendLine(R"(plot "data2.dat" with lines)");
    gp.sendLine(R"(plot "conv2.dat" with lines)");
    gp.sendLine(R"(plot "conv8.dat" with lines)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "conv3.dat" with lines)");
    gp.sendLine(R"(plot "data13.dat" with lines)");

    //gp.sendLine(R"(plot "data1.dat" with lines)");
    //gp.sendLine(R"(plot "data9.dat" with lines)");
    //gp.sendLine(R"(plot "conv5.dat" with lines)");
    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "conv4.dat" with lines)");
    gp.sendLine(R"(plot "data12.dat" with lines)");

    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "data10.dat" with lines)");
    gp.sendLine(R"(plot "conv6.dat" with lines)");
    gp.sendLine(R"(plot "conv10.dat" with lines)");

    return 0;
}