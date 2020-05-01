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

    ComplexPlot frequencyStep_cutted;              // step
    ComplexPlot frequencyStep_transform_cutted;    // step -> frequency spectra

    ComplexPlot sum_sins;
    ComplexPlot sum_sins_FILTERED;
    ComplexPlot sum_sins_SPECTRA1;          // sin(x) -> frequency spectra
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
    double filter_width = (double)quantity/4;

    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;
        double t = (double)i/quantity * 10 * PI;

        double gaussWindow = 0;
        if(!(i > quantity/2 - filter_width && i < quantity/2 + filter_width))
            gaussWindow = exp(-(i-quantity/2 + filter_width)*(i-quantity/2 + filter_width)/100) + exp(-(i-quantity/2-filter_width)*(i-quantity/2-filter_width)/100);
        // quantity - общее число точек
        // filter_width - ширина фильтра (измеряемая в количестве точек)
        double phase = 0; // По умолчанию фаза = 0
        if(i < filter_width)
            // Если мы в начале, то фаза лежит в интервале [0, PI]
            phase = 0*PI*i/filter_width;
        else if(i > quantity - filter_width)
            // Если мы в начале, то фаза лежит в интервале [-PI, 0]
            phase = 0*PI*(-1 + (i - quantity + filter_width)/filter_width);

        if(i > quantity/2 - filter_width && i < quantity/2 + filter_width)
            // Если мы в начале или в конце, то задаём модуль единицу, учитывая фазу:
            frequencyStep.push(i, {1.0001*cos(phase), 1.0001*sin(phase)});
        else
            // Если мы в середине, то модуль ноль и фаза ноль.
            frequencyStep.push(i, {gaussWindow*cos(phase), gaussWindow*sin(phase)});

        sum_sins.push(i, {sin(30*t) + sin(4*t), 0});

        sin_high.push(i, {sin(30*t), 0});

        sin_low.push(i, {sin(x), 0});

        if(i < quantity/2)
            half_sin_low_half_sin_high.push(i, {sin(4*t), 0});
        else
            half_sin_low_half_sin_high.push(i, {sin(30*t), 0});

        if((int)(i / 10) % 2 == 0)
            periodic_step.push(i, {1.0001, 0});
        else
            periodic_step.push(i, {0, 0});
    }

    int Ncut = quantity/20;
    fourierTransform(frequencyStep, frequencyStep_transform_cutted);
    frequencyStep_transform_cutted.cut(Ncut, frequencyStep_transform_cutted.size() - Ncut );
    inverseFourierTransform(frequencyStep_transform_cutted, frequencyStep_cutted);

    addNoise(half_sin_low_half_sin_high);
    addNoise(sum_sins, 0.3);

    // SIMPLE METHOD
    ComplexPlot resCross;
    ComplexPlot inv_resCross;

    ComplexPlot conv;
    ComplexPlot conv_inv;

    //fftw_fourierTransform(frequencyStep, frequencyStep_transform);
    fourierTransform(frequencyStep, frequencyStep_transform);
    //inverseFourierTransform(frequencyStep, frequencyStep_transform);
    // sum sins
    fourierTransform(sum_sins, sum_sins_SPECTRA1);
    ComplexPlot cross_sum_sins;
    //convolution(frequencyStep_transform, sum_sins, sum_sins_FILTERED);
    cross(sum_sins_SPECTRA1, frequencyStep_cutted, cross_sum_sins);
    inverseFourierTransform(cross_sum_sins, sum_sins_FILTERED);
    fourierTransform(sum_sins_FILTERED, sum_sins_SPECTRA2);

    // sin high
    fourierTransform(sin_high, sin_high_SPECTRA1);
    ComplexPlot cross_sin_high;
    //convolution(frequencyStep_transform, sin_high, sin_high_FILTERED);
    cross(sin_high_SPECTRA1, frequencyStep_cutted, cross_sin_high);
    inverseFourierTransform(cross_sin_high, sin_high_FILTERED);
    fourierTransform(sin_high_FILTERED, sin_high_SPECTRA2);

    // sin low
    fourierTransform(sin_low, sin_low_SPECTRA1);
    ComplexPlot cross_sin_low;
    //convolution(frequencyStep_transform, sin_low, sin_low_FILTERED);
    cross(sin_low_SPECTRA1, frequencyStep_cutted, cross_sin_low);
    inverseFourierTransform(cross_sin_low, sin_low_FILTERED);
    fourierTransform(sin_low_FILTERED, sin_low_SPECTRA2);

    // half sin low, half sin high
    fourierTransform(half_sin_low_half_sin_high, half_sin_low_half_sin_high_SPECTRA1);
    ComplexPlot cross_half_sin_low_half_sin_high;
    //convolution(frequencyStep_transform, half_sin_low_half_sin_high, half_sin_low_half_sin_high_FILTERED);
    cross(half_sin_low_half_sin_high_SPECTRA1, frequencyStep_cutted, cross_half_sin_low_half_sin_high);
    inverseFourierTransform(cross_half_sin_low_half_sin_high, half_sin_low_half_sin_high_FILTERED);
    fourierTransform(half_sin_low_half_sin_high_FILTERED, half_sin_low_half_sin_high_SPECTRA2);

    // periodic step
    fourierTransform(periodic_step, periodic_step_SPECTRA1);
    ComplexPlot cross_periodic_step;
    //convolution(frequencyStep_transform, periodic_step, periodic_step_FILTERED);
    cross(periodic_step_SPECTRA1, frequencyStep_cutted, cross_periodic_step);
    inverseFourierTransform(cross_periodic_step, periodic_step_FILTERED);
    fourierTransform(periodic_step_FILTERED, periodic_step_SPECTRA2);

    // SAVING DATA

    saveVectorPoint2DToFile(frequencyStep.abs(), "frequencyStep.dat");
    saveVectorPoint2DToFile(frequencyStep.phase(), "frequencyStep_phase.dat");

    saveVectorPoint2DToFile(frequencyStep_cutted.abs(), "frequencyStep_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_cutted.phase(), "frequencyStep_cutted_phase.dat");

    saveVectorPoint2DToFile(frequencyStep_transform.abs(), "frequencyStep_transform.dat");
    saveVectorPoint2DToFile(frequencyStep_transform.phase(), "frequencyStep_transform_phase.dat");

    saveVectorPoint2DToFile(frequencyStep_transform_cutted.abs(), "frequencyStep_transform_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_transform_cutted.phase(), "frequencyStep_transform_cutted_phase.dat");

    saveVectorPoint2DToFile(frequencyStep_transform.real(), "frequencyStep_transform_real.dat");
    saveVectorPoint2DToFile(frequencyStep_transform.imagine(), "frequencyStep_transform_image.dat");

    saveVectorPoint2DToFile(frequencyStep_transform_cutted.real(), "frequencyStep_transform_cutted_real.dat");
    saveVectorPoint2DToFile(frequencyStep_transform_cutted.imagine(), "frequencyStep_transform_cutted_image.dat");

    saveVectorPoint2DToFile(sum_sins.real(), "sum_sins.dat");
    saveVectorPoint2DToFile(sum_sins_FILTERED.real(), "sum_sins_FILTERED.dat", sum_sins_FILTERED.size());
    saveVectorPoint2DToFile(sum_sins_SPECTRA1.abs(), "sum_sins_SPECTRA1.dat");
    saveVectorPoint2DToFile(sum_sins_SPECTRA2.abs(), "sum_sins_SPECTRA2.dat");

    saveVectorPoint2DToFile(sin_high.real(), "sin_high.dat");
    saveVectorPoint2DToFile(sin_high_FILTERED.real(), "sin_high_FILTERED.dat", sin_high_FILTERED.size());
    saveVectorPoint2DToFile(sin_high_SPECTRA1.abs(), "sin_high_SPECTRA1.dat");
    saveVectorPoint2DToFile(sin_high_SPECTRA2.abs(), "sin_high_SPECTRA2.dat");

    saveVectorPoint2DToFile(sin_low.real(), "sin_low.dat");
    saveVectorPoint2DToFile(sin_low_FILTERED.real(), "sin_low_FILTERED.dat", sin_low_FILTERED.size());
    saveVectorPoint2DToFile(sin_low_SPECTRA1.abs(), "sin_low_SPECTRA1.dat");
    saveVectorPoint2DToFile(sin_low_SPECTRA2.abs(), "sin_low_SPECTRA2.dat");

    saveVectorPoint2DToFile(half_sin_low_half_sin_high.real(), "half_sin_low_half_sin_high.dat");
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_FILTERED.real(), "half_sin_low_half_sin_high_FILTERED.dat", half_sin_low_half_sin_high_FILTERED.size());
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_SPECTRA1.abs(), "half_sin_low_half_sin_high_SPECTRA1.dat");
    saveVectorPoint2DToFile(half_sin_low_half_sin_high_SPECTRA2.abs(), "half_sin_low_half_sin_high_SPECTRA2.dat");

    saveVectorPoint2DToFile(periodic_step.real(), "periodic_step.dat");
    saveVectorPoint2DToFile(periodic_step_FILTERED.real(), "periodic_step_FILTERED.dat", periodic_step_FILTERED.size());
    saveVectorPoint2DToFile(periodic_step_SPECTRA1.abs(), "periodic_step_SPECTRA1.dat");
    saveVectorPoint2DToFile(periodic_step_SPECTRA2.abs(), "periodic_step_SPECTRA2.dat");

    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 4, 1)");
    gp.sendLine(R"(set key spacing 1.5)");
    gp.sendLine(R"(unset key)");

    // step check:
    //gp.sendLine(R"(plot "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_cutted_phase.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_transform_cutted_real.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_transform_cutted_image.dat" with lines)");


    //gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_SPECTRA2.dat" with lines)");

    //gp.sendLine(R"(plot "sin_high.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "sin_low.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA2.dat" with lines)");
//
    gp.sendLine(R"(plot "periodic_step.dat" with lines)");
    gp.sendLine(R"(plot "periodic_step_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    gp.sendLine(R"(plot "periodic_step_FILTERED.dat" with lines)");
    gp.sendLine(R"(plot "periodic_step_SPECTRA2.dat" with lines)");

    return 0;
}