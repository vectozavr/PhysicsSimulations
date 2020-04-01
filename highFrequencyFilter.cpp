//
// Created by ivan- on 15.03.2020.
//

#include <algorithm>
#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;
using namespace std;

int main() {
    ComplexPlot frequencyStep;              // step
    ComplexPlot frequencyStep_transform;    // step -> frequency spectra

    ComplexPlot data1;          // sin(x)
    ComplexPlot transform1;     // sin(x) -> frequency spectra


    int quantity = 500;
    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;

        if(i < quantity/10)
            frequencyStep.push(i, {1, 0});
        else
            frequencyStep.push(i, {0, 0});

        data1.push(i, {sin(x) + sin(4*x)/5, 0});
    }

    addNoise(data1);

    // SIMPLE METHOD
    ComplexPlot resCross;
    ComplexPlot inv_resCross;

    fourierTransform(frequencyStep, frequencyStep_transform);
    fourierTransform(data1, transform1);
    cross(frequencyStep, transform1, resCross);
    inverseFourierTransform(resCross, inv_resCross);
    // 2 METHOD
    ComplexPlot resConv;
    ComplexPlot test;

    fourierTransform(frequencyStep, frequencyStep_transform);
    convolution(data1, frequencyStep_transform, resConv);
    convolution(frequencyStep, frequencyStep, test);
    // SAVING DATA
    saveVectorPoint2DToFile(data1.real(), "data1.dat");
    saveVectorPoint2DToFile(frequencyStep.real(), "frequencyStep.dat");
    saveVectorPoint2DToFile(transform1.abs(), "transform1.dat");
    saveVectorPoint2DToFile(inv_resCross.real(), "inv_resCross.dat");

    saveVectorPoint2DToFile(frequencyStep_transform.real(), "frequencyStep_transform.dat");
    resConv.cut(resConv.size()/2, resConv.size());
    saveVectorPoint2DToFile(resConv.real(), "resConv.dat", resConv.size()/2);

    saveVectorPoint2DToFile(test.real(), "test.dat");
    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 2, 3)");
    gp.sendLine(R"(unset key)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "frequencyStep.dat" with lines)");
    gp.sendLine(R"(plot "inv_resCross.dat" with lines)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "frequencyStep_transform.dat" with lines)");
    gp.sendLine(R"(plot "resConv.dat" with lines)");

    return 0;
}