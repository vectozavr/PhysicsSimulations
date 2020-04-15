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

    ComplexPlot sum_sins;
    ComplexPlot sum_sins_FILTERED;
    ComplexPlot sum_sins_SPECTRA1;     // sin(x) -> frequency spectra
    ComplexPlot sum_sins_SPECTRA2;

    ComplexPlot sin_high;
    ComplexPlot sin_low;
    ComplexPlot half_sin_low_half_sin_high;
    ComplexPlot periodic_step;

    ComplexPlot sin_high_FILTERED;
    ComplexPlot sin_low_FILTERED;
    ComplexPlot half_sin_low_half_sin_high_FILTERED;
    ComplexPlot periodic_step_FILTERED;

    ComplexPlot sin_high_SPECTRA1;
    ComplexPlot sin_low_SPECTRA1;
    ComplexPlot half_sin_low_half_sin_high_SPECTRA1;
    ComplexPlot periodic_step_SPECTRA1;

    ComplexPlot sin_high_SPECTRA2;
    ComplexPlot sin_low_SPECTRA2;
    ComplexPlot half_sin_low_half_sin_high_SPECTRA2;
    ComplexPlot periodic_step_SPECTRA2;

    int quantity = 500;
    double filter_width = (double)quantity/10;

    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;
        double t = (double)i/quantity * 10 * PI;

        double phase = 0;
        if(i < filter_width)
            phase = 2*PI*i/filter_width;
        else if(i > quantity - filter_width)
            phase = 2*PI*(-1 + (i - quantity + filter_width)/filter_width);


        if(i < filter_width || i > quantity - filter_width)
            frequencyStep.push(i, {1*cos(phase), 1*sin(phase)});
        else
            frequencyStep.push(i, {0, 0});

        sum_sins.push(i, {sin(t) + sin(4*t), 0});

        sin_high.push(i, {sin(15*t), 0});

        sin_low.push(i, {sin(x), 0});

        if(i < quantity/2)
            half_sin_low_half_sin_high.push(i, {sin(x), 0});
        else
            half_sin_low_half_sin_high.push(i, {sin(10*x), 0});

        if((int)(i / 10) % 2 == 0)
            periodic_step.push(i, {1.0001, 0});
        else
            periodic_step.push(i, {0, 0});
    }

    //addNoise(sum_sins);

    // SIMPLE METHOD
    ComplexPlot resCross;
    ComplexPlot inv_resCross;

    ComplexPlot conv;
    ComplexPlot conv_inv;

    inverseFourierTransform(frequencyStep, frequencyStep_transform);

    convolution(frequencyStep_transform, sum_sins, sum_sins_FILTERED);
    fourierTransform(sum_sins_FILTERED, sum_sins_SPECTRA2);
    fourierTransform(sum_sins, sum_sins_SPECTRA1);

    convolution(frequencyStep_transform, sin_high, sin_high_FILTERED);
    fourierTransform(sin_high_FILTERED, sin_high_SPECTRA2);
    fourierTransform(sin_high, sin_high_SPECTRA1);

    convolution(frequencyStep_transform, sin_low, sin_low_FILTERED);
    fourierTransform(sin_low_FILTERED, sin_low_SPECTRA2);
    fourierTransform(sin_low, sin_low_SPECTRA1);
//
    convolution(frequencyStep_transform, half_sin_low_half_sin_high, half_sin_low_half_sin_high_FILTERED);
    fourierTransform(half_sin_low_half_sin_high_FILTERED, half_sin_low_half_sin_high_SPECTRA2);
    fourierTransform(half_sin_low_half_sin_high, half_sin_low_half_sin_high_SPECTRA1);
//
    convolution(frequencyStep_transform, periodic_step, periodic_step_FILTERED);
    fourierTransform(periodic_step_FILTERED, periodic_step_SPECTRA2);
    fourierTransform(periodic_step, periodic_step_SPECTRA1);

    // SAVING DATA

    saveVectorPoint2DToFile(frequencyStep.abs(), "frequencyStep.dat");
    saveVectorPoint2DToFile(frequencyStep.phase(), "frequencyStep_phase.dat");

    saveVectorPoint2DToFile(frequencyStep_transform.abs(), "frequencyStep_transform.dat");
    saveVectorPoint2DToFile(frequencyStep_transform.phase(), "frequencyStep_transform_phase.dat");

    saveVectorPoint2DToFile(sum_sins.real(), "sum_sins.dat");
    saveVectorPoint2DToFile(sum_sins_FILTERED.real(), "sum_sins_FILTERED.dat", sum_sins_FILTERED.size()/2);
    saveVectorPoint2DToFile(sum_sins_SPECTRA1.abs(), "sum_sins_SPECTRA1.dat");
    saveVectorPoint2DToFile(sum_sins_SPECTRA2.abs(), "sum_sins_SPECTRA2.dat");

    saveVectorPoint2DToFile(sin_high.real(), "sin_high.dat");
    saveVectorPoint2DToFile(sin_high_FILTERED.real(), "sin_high_FILTERED.dat", sin_high_FILTERED.size()/2);
    saveVectorPoint2DToFile(sin_high_SPECTRA1.abs(), "sin_high_SPECTRA1.dat");
    saveVectorPoint2DToFile(sin_high_SPECTRA2.abs(), "sin_high_SPECTRA2.dat");

    saveVectorPoint2DToFile(sin_low.real(), "sin_low.dat");
    saveVectorPoint2DToFile(sin_low_FILTERED.real(), "sin_low_FILTERED.dat", sin_low_FILTERED.size()/2);
    saveVectorPoint2DToFile(sin_low_SPECTRA1.abs(), "sin_low_SPECTRA1.dat");
    saveVectorPoint2DToFile(sin_low_SPECTRA2.abs(), "sin_low_SPECTRA2.dat");

    saveVectorPoint2DToFile(half_sin_low_half_sin_high.real(), "half_sin_low_half_sin_high.dat");
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_FILTERED.real(), "half_sin_low_half_sin_high_FILTERED.dat", half_sin_low_half_sin_high_FILTERED.size()/2);
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_SPECTRA1.abs(), "half_sin_low_half_sin_high_SPECTRA1.dat");
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_SPECTRA2.abs(), "half_sin_low_half_sin_high_SPECTRA2.dat");

    saveVectorPoint2DToFile(periodic_step.real(), "periodic_step.dat");
    saveVectorPoint2DToFile(periodic_step_FILTERED.real(), "periodic_step_FILTERED.dat", periodic_step_FILTERED.size()/2);
    saveVectorPoint2DToFile(periodic_step_SPECTRA1.abs(), "periodic_step_SPECTRA1.dat");
    saveVectorPoint2DToFile(periodic_step_SPECTRA2.abs(), "periodic_step_SPECTRA2.dat");

    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 4, 1)");
    gp.sendLine(R"(set key spacing 1.5)");
    gp.sendLine(R"(unset key)");


    //gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_phase.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_SPECTRA2.dat" with lines)");

    gp.sendLine(R"(plot "sin_high.dat" with lines)");
    gp.sendLine(R"(plot "sin_high_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    gp.sendLine(R"(plot "sin_high_FILTERED.dat" with lines)");
    gp.sendLine(R"(plot "sin_high_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "sin_low.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "periodic_step.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_SPECTRA2.dat" with lines)");

    return 0;
}