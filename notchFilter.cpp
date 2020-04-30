//
// Created by ivan- on 29.04.2020.
//

#include <algorithm>
#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;
using namespace std;

int main() {

    ComplexPlot frequencyStep_high;              // step
    ComplexPlot frequencyStep_high_transform;    // step -> frequency spectra
    ComplexPlot frequencyStep_high_cutted;              // step
    ComplexPlot frequencyStep_high_transform_cutted;    // step -> frequency spectra

    ComplexPlot frequencyStep_low;              // step
    ComplexPlot frequencyStep_low_transform;    // step -> frequency spectra
    ComplexPlot frequencyStep_low_cutted;              // step
    ComplexPlot frequencyStep_low_transform_cutted;    // step -> frequency spectra

    ComplexPlot frequencyStep_high_low_cutted;

    ComplexPlot sum_sins;
    ComplexPlot sum_sins_FILTERED;
    ComplexPlot sum_sins_FILTERED_REG;
    ComplexPlot sum_sins_SPECTRA1;          // sin(x) -> frequency spectra
    ComplexPlot sum_sins_SPECTRA2;

    ComplexPlot sum_sins_FILTERED_REG_SPECTRA;

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
    double filter_high_width = (double)quantity/6;
    double filter_low_width = (double)quantity/10;

    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;
        double t = (double)i/quantity * 10 * PI;

        double gaussWindow_high = 0;
        double gaussWindow_low = 0;
        if(!(i > quantity/2 - gaussWindow_high && i < quantity/2 + gaussWindow_high))
            gaussWindow_high = exp(-(i-quantity/2 + gaussWindow_high)*(i-quantity/2 + gaussWindow_high)/100) + exp(-(i-quantity/2-gaussWindow_high)*(i-quantity/2-gaussWindow_high)/100);
        if((i >= gaussWindow_low) && (i <= quantity - gaussWindow_low))
            gaussWindow_low =  exp(-(i-gaussWindow_low)*(i-gaussWindow_low)/100) + exp(-(i-quantity+gaussWindow_low)*(i-quantity+gaussWindow_low)/100);
        // quantity - общее число точек
        // filter_width - ширина фильтра (измеряемая в количестве точек)
        double phase = 0; // По умолчанию фаза = 0

        if(i > quantity/2 - filter_high_width && i < quantity/2 + filter_high_width)
            frequencyStep_high.push(i, {1.0001*cos(phase), 1.0001*sin(phase)});
        else
            frequencyStep_high.push(i, {gaussWindow_high*cos(phase), gaussWindow_high*sin(phase)});

        if(i < filter_low_width || i > quantity - filter_low_width)
            frequencyStep_low.push(i, {1.0001*cos(phase), 1.0001*sin(phase)});
        else
            frequencyStep_low.push(i, {gaussWindow_low*cos(phase), gaussWindow_low*sin(phase)});

        sum_sins.push(i, {sin(30*t) + sin(4*t) + sin(8*t) + sin(16*t), 0});

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
    fourierTransform(frequencyStep_high, frequencyStep_high_transform_cutted);
    frequencyStep_high_transform_cutted.cut(Ncut, frequencyStep_high_transform_cutted.size() - Ncut );
    inverseFourierTransform(frequencyStep_high_transform_cutted, frequencyStep_high_cutted);

    fourierTransform(frequencyStep_low, frequencyStep_low_transform_cutted);
    frequencyStep_low_transform_cutted.cut(Ncut, frequencyStep_low_transform_cutted.size() - Ncut );
    inverseFourierTransform(frequencyStep_low_transform_cutted, frequencyStep_low_cutted);

    for(int k = 0; k < frequencyStep_high_cutted.size(); k++)
        frequencyStep_high_low_cutted.push(k, frequencyStep_high_cutted.v_c[k].second + frequencyStep_low_cutted.v_c[k].second);

    //addNoise(half_sin_low_half_sin_high);
    //addNoise(sum_sins, 0.3);

    // SIMPLE METHOD
    ComplexPlot resCross;
    ComplexPlot inv_resCross;

    ComplexPlot conv;
    ComplexPlot conv_inv;

    //fftw_fourierTransform(frequencyStep, frequencyStep_transform);
    fourierTransform(frequencyStep_high, frequencyStep_high_transform);
    //inverseFourierTransform(frequencyStep, frequencyStep_transform);

    // sum sins
    // high frequency filter
    fourierTransform(sum_sins, sum_sins_SPECTRA1);
    ComplexPlot cross_sum_sins;
    //convolution(frequencyStep_transform, sum_sins, sum_sins_FILTERED);
    cross(sum_sins_SPECTRA1, frequencyStep_high_low_cutted, cross_sum_sins);
    inverseFourierTransform(cross_sum_sins, sum_sins_FILTERED);
    fourierTransform(sum_sins_FILTERED, sum_sins_SPECTRA2);

    for(int i = 0; i < sum_sins_FILTERED.size(); i++)
        sum_sins_FILTERED_REG.push(i, sum_sins.v_c[i].second - sum_sins_FILTERED.v_c[i].second);

    fourierTransform(sum_sins_FILTERED_REG, sum_sins_FILTERED_REG_SPECTRA);

    // low frequency filter

    // sin high
    // high frequency filter
    fourierTransform(sin_high, sin_high_SPECTRA1);
    ComplexPlot cross_sin_high;
    //convolution(frequencyStep_transform, sin_high, sin_high_FILTERED);
    cross(sin_high_SPECTRA1, frequencyStep_high_cutted, cross_sin_high);
    inverseFourierTransform(cross_sin_high, sin_high_FILTERED);
    fourierTransform(sin_high_FILTERED, sin_high_SPECTRA2);

    // low frequency filter

    // sin low
    fourierTransform(sin_low, sin_low_SPECTRA1);
    ComplexPlot cross_sin_low;
    //convolution(frequencyStep_transform, sin_low, sin_low_FILTERED);
    cross(sin_low_SPECTRA1, frequencyStep_high_cutted, cross_sin_low);
    inverseFourierTransform(cross_sin_low, sin_low_FILTERED);
    fourierTransform(sin_low_FILTERED, sin_low_SPECTRA2);

    // half sin low, half sin high
    fourierTransform(half_sin_low_half_sin_high, half_sin_low_half_sin_high_SPECTRA1);
    ComplexPlot cross_half_sin_low_half_sin_high;
    //convolution(frequencyStep_transform, half_sin_low_half_sin_high, half_sin_low_half_sin_high_FILTERED);
    cross(half_sin_low_half_sin_high_SPECTRA1, frequencyStep_high_cutted, cross_half_sin_low_half_sin_high);
    inverseFourierTransform(cross_half_sin_low_half_sin_high, half_sin_low_half_sin_high_FILTERED);
    fourierTransform(half_sin_low_half_sin_high_FILTERED, half_sin_low_half_sin_high_SPECTRA2);

    // periodic step
    fourierTransform(periodic_step, periodic_step_SPECTRA1);
    ComplexPlot cross_periodic_step;
    //convolution(frequencyStep_transform, periodic_step, periodic_step_FILTERED);
    cross(periodic_step_SPECTRA1, frequencyStep_high_cutted, cross_periodic_step);
    inverseFourierTransform(cross_periodic_step, periodic_step_FILTERED);
    fourierTransform(periodic_step_FILTERED, periodic_step_SPECTRA2);

    // SAVING DATA
    saveVectorPoint2DToFile(frequencyStep_high.abs(), "frequencyStep_high.dat");
    saveVectorPoint2DToFile(frequencyStep_high.phase(), "frequencyStep_high_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_high_cutted.abs(), "frequencyStep_high_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_high_cutted.phase(), "frequencyStep_high_cutted_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform.abs(), "frequencyStep_high_transform.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform.phase(), "frequencyStep_high_transform_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform_cutted.abs(), "frequencyStep_high_transform_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform_cutted.phase(), "frequencyStep_high_transform_cutted_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform.real(), "frequencyStep_high_transform_real.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform.imagine(), "frequencyStep_high_transform_image.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform_cutted.real(), "frequencyStep_high_transform_cutted_real.dat");
    saveVectorPoint2DToFile(frequencyStep_high_transform_cutted.imagine(), "frequencyStep_high_transform_cutted_image.dat");

    saveVectorPoint2DToFile(frequencyStep_low.abs(), "frequencyStep_low.dat");
    saveVectorPoint2DToFile(frequencyStep_low.phase(), "frequencyStep_low_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_low_cutted.abs(), "frequencyStep_low_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_low_cutted.phase(), "frequencyStep_low_cutted_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform.abs(), "frequencyStep_low_transform.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform.phase(), "frequencyStep_low_transform_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform_cutted.abs(), "frequencyStep_low_transform_cutted.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform_cutted.phase(), "frequencyStep_low_transform_cutted_phase.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform.real(), "frequencyStep_low_transform_real.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform.imagine(), "frequencyStep_low_transform_image.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform_cutted.real(), "frequencyStep_low_transform_cutted_real.dat");
    saveVectorPoint2DToFile(frequencyStep_low_transform_cutted.imagine(), "frequencyStep_low_transform_cutted_image.dat");


    saveVectorPoint2DToFile(sum_sins.real(), "sum_sins.dat");
    saveVectorPoint2DToFile(sum_sins_FILTERED_REG.real(), "sum_sins_FILTERED_REG.dat");
    saveVectorPoint2DToFile(sum_sins_FILTERED_REG_SPECTRA.abs(), "sum_sins_FILTERED_REG_SPECTRA.dat");
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
    //gp.sendLine(R"(plot "frequencyStep_high_cutted.dat" with lines,                 "frequencyStep_low_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_high_cutted_phase.dat" with lines,           "frequencyStep_low_cutted_phase.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_high_transform_cutted_real.dat" with lines,  "frequencyStep_low_transform_cutted_real.dat" with lines)");
    //gp.sendLine(R"(plot "frequencyStep_high_transform_cutted_image.dat" with lines, "frequencyStep_low_transform_cutted_image.dat" with lines)");

    // УЗКОПОЛОСНЫЙ ФИЛЬТР РЕЗУЛЬТАТ
    gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_SPECTRA1.dat" with lines, "frequencyStep_high_cutted.dat" with lines, "frequencyStep_low_cutted.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_FILTERED.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_SPECTRA2.dat" with lines)");
    // РЕЖЕКТОРНЫЙ ФИЛЬТР РЕЗУЛЬТАТ
    //gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_SPECTRA1.dat" with lines, "frequencyStep_high_cutted.dat" with lines, "frequencyStep_low_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_FILTERED_REG.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_FILTERED_REG_SPECTRA.dat" with lines)");

    //gp.sendLine(R"(plot "sin_high.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sin_high_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "sin_low.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "sin_low_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA1.dat" with lines, "frequencyStep_cutted.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "half_sin_low_half_sin_high_SPECTRA2.dat" with lines)");
//
    //gp.sendLine(R"(plot "periodic_step.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_SPECTRA1.dat" with lines, "frequencyStep.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_FILTERED.dat" with lines)");
    //gp.sendLine(R"(plot "periodic_step_SPECTRA2.dat" with lines)");

    return 0;
}