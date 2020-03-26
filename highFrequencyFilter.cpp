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

        if(i < quantity/2)
            frequencyStep.push(i, {1, 0});
        else
            frequencyStep.push(i, {0, 0});

        data1.push(i, {sin(x), 0});
    }

    ComplexPlot resultConv1;

    fourierTransform(frequencyStep, frequencyStep_transform);
    fourierTransform(data1, transform1);
    convolution(data1, frequencyStep, resultConv1);
    //inverseFourierTransform(resultConv1, resultConv_inv1);

    saveVectorPoint2DToFile(data1.real(), "data1.dat");
    saveVectorPoint2DToFile(transform1.abs(), "transform1.dat");
    saveVectorPoint2DToFile(resultConv1.real(), "resultConv1.dat");
    saveVectorPoint2DToFile(frequencyStep.real(), "frequencyStep.dat");
    saveVectorPoint2DToFile(resultConv1.real(), "resultConv1.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");
    //gp.sendLine(R"(set multiplot layout 2, 1)");
    gp.sendLine(R"(unset key)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep.dat" with lines)");
    return 0;
}