//
// Created by ivan- on 01.05.2020.
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

    int quantity = 500;
    double filter_width = (double)quantity/10;

    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;
        double t = (double)i/quantity * 10 * PI;

        double gaussWindow = 0;
        if((i >= filter_width) && (i <= quantity - filter_width))
            gaussWindow = exp(-(i-filter_width)*(i-filter_width)/100) + exp(-(i-quantity+filter_width)*(i-quantity+filter_width)/100);
        // quantity - общее число точек
        // filter_width - ширина фильтра (измеряемая в количестве точек)
        double phase = 0; // По умолчанию фаза = 0
        if(i < filter_width)
            // Если мы в начале, то фаза лежит в интервале [0, PI]
            phase = 0*PI*i/filter_width;
        else if(i > quantity - filter_width)
            // Если мы в начале, то фаза лежит в интервале [-PI, 0]
            phase = 0*PI*(-1 + (i - quantity + filter_width)/filter_width);

        if(i < filter_width || i > quantity - filter_width)
            // Если мы в начале или в конце, то задаём модуль единицу, учитывая фазу:
            frequencyStep.push(i, {1.0001*cos(phase), 1.0001*sin(phase)});
        else
            // Если мы в середине, то модуль ноль и фаза ноль.
            frequencyStep.push(i, {gaussWindow*cos(phase), gaussWindow*sin(phase)});

        sum_sins.push(i, {sin(t) + sin(4*t), 0});
    }

    int P = 4;
    int Q = 5;

    // Increase frequency in P times
    ComplexPlot sum_sins_P;
    ComplexPlot sum_sins_P_SPECTRA;

    ComplexPlot sum_sins_P_FILTERED;
    ComplexPlot sum_sins_P_FILTERED_SPECTRA;

    ComplexPlot freq_Step_1;
    ComplexPlot freq_Step_1_SPECTRA;
    ComplexPlot freq_Step_1_CUTTED;
    int indexP = 0;
    for(int i = 0; i < sum_sins.size(); i++) {
        sum_sins_P.push(indexP++, sum_sins.v_c[i].second);
        for(int k = 0; k < P-1; k++) {
            sum_sins_P.push(indexP++, 0);
        }
    }
    fourierTransform(sum_sins_P, sum_sins_P_SPECTRA);
    for(int i = 0; i < indexP; i++) {
        double gaussWindow = 0;
        if((i >= (double)sum_sins.size()/2) && (i <= indexP - (double)sum_sins.size()/2))
            gaussWindow = exp(-(i-(double)sum_sins.size()/2)*(i-(double)sum_sins.size()/2)/100) + exp(-(i-indexP+(double)sum_sins.size()/2)*(i-indexP+(double)sum_sins.size()/2)/100);;

        double phase = 0; // По умолчанию фаза = 0

        if(i < (double)sum_sins.size()/2 || i > indexP - (double)sum_sins.size()/2)
            freq_Step_1.push(i, {1.0001*P*cos(phase), P*sin(phase)});
        else
            freq_Step_1.push(i, {1.0001*P*gaussWindow*cos(phase), P*gaussWindow*sin(phase)});
    }

    int Ncut1 = indexP/20;
    fourierTransform(freq_Step_1, freq_Step_1_SPECTRA);
    freq_Step_1_SPECTRA.cut(Ncut1, freq_Step_1_SPECTRA.size() - Ncut1 );
    inverseFourierTransform(freq_Step_1_SPECTRA, freq_Step_1_CUTTED);

    ComplexPlot cross_sum_sins;
    cross(sum_sins_P_SPECTRA, freq_Step_1_CUTTED, cross_sum_sins);
    inverseFourierTransform(cross_sum_sins, sum_sins_P_FILTERED);
    fourierTransform(sum_sins_P_FILTERED, sum_sins_P_FILTERED_SPECTRA);

    // Decreasing frequency in Q times
    ComplexPlot sum_sins_P_Q_RESULT;
    ComplexPlot sum_sins_P_Q_SPECTRA;

    ComplexPlot sum_sins_P_Q_FILTERED;
    ComplexPlot sum_sins_P_Q_FILTERED_SPECTRA;

    ComplexPlot freq_Step_2;
    ComplexPlot freq_Step_2_SPECTRA;
    ComplexPlot freq_Step_2_CUTTED;



    for(int i = 0; i < (double)sum_sins_P_FILTERED.size(); i++) {
        double gaussWindow = 0;
        if((i >= (double)sum_sins_P_FILTERED.size()/(2*Q)) && (i <= (double)sum_sins_P_FILTERED.size() - (double)sum_sins_P_FILTERED.size()/(2*Q)))
            gaussWindow = exp(-(i-(double)sum_sins_P_FILTERED.size()/(2*Q))*(i-(double)sum_sins_P_FILTERED.size()/(2*Q))/100) + exp(-(i-(double)sum_sins_P_FILTERED.size()+(double)sum_sins_P_FILTERED.size()/(2*Q))*(i-(double)sum_sins_P_FILTERED.size()+(double)sum_sins_P_FILTERED.size()/(2*Q))/100);;

        double phase = 0; // По умолчанию фаза = 0

        if(i < (double)sum_sins_P_FILTERED.size()/(2*Q) || i > (double)sum_sins_P_FILTERED.size() - (double)sum_sins_P_FILTERED.size()/(2*Q))
            freq_Step_2.push(i, {1.0001*cos(phase), 1.0001*sin(phase)});
        else
            freq_Step_2.push(i, {1.0001*gaussWindow*cos(phase), 1.0001*gaussWindow*sin(phase)});
    }

    int Ncut2 = (double)freq_Step_2.size()/20;
    fourierTransform(freq_Step_2, freq_Step_2_SPECTRA);
    freq_Step_2_SPECTRA.cut(Ncut2, freq_Step_2_SPECTRA.size() - Ncut2 );
    inverseFourierTransform(freq_Step_2_SPECTRA, freq_Step_2_CUTTED);

    ComplexPlot cross_sum_sins2;
    cross(sum_sins_P_FILTERED_SPECTRA, freq_Step_2_CUTTED, cross_sum_sins2);
    inverseFourierTransform(cross_sum_sins2, sum_sins_P_Q_FILTERED);
    fourierTransform(sum_sins_P_Q_FILTERED, sum_sins_P_Q_FILTERED_SPECTRA);

    int indexQ = 0;
    for(int i = 0; i < sum_sins_P_Q_FILTERED.size()/Q; i++) {
        sum_sins_P_Q_RESULT.push(indexQ++, sum_sins_P_FILTERED.v_c[i*Q].second);
    }
    fourierTransform(sum_sins_P_Q_RESULT, sum_sins_P_Q_SPECTRA);

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

    saveVectorPoint2DToFile(sum_sins_P.real(), "sum_sins_P.dat");
    saveVectorPoint2DToFile(sum_sins_P_SPECTRA.abs(), "sum_sins_P_SPECTRA.dat");
    saveVectorPoint2DToFile(freq_Step_1.abs(), "freq_Step_1.dat");
    saveVectorPoint2DToFile(freq_Step_1_CUTTED.abs(), "freq_Step_1_CUTTED.dat");
    saveVectorPoint2DToFile(sum_sins_P_FILTERED.real(), "sum_sins_P_FILTERED.dat");
    saveVectorPoint2DToFile(sum_sins_P_FILTERED_SPECTRA.abs(), "sum_sins_P_FILTERED_SPECTRA.dat");

    saveVectorPoint2DToFile(sum_sins_P_Q_RESULT.real(), "sum_sins_P_Q_RESULT.dat");
    saveVectorPoint2DToFile(sum_sins_P_Q_SPECTRA.abs(), "sum_sins_P_Q_SPECTRA.dat");
    saveVectorPoint2DToFile(freq_Step_2.abs(), "freq_Step_2.dat");
    saveVectorPoint2DToFile(freq_Step_2_CUTTED.abs(), "freq_Step_2_CUTTED.dat");
    saveVectorPoint2DToFile(sum_sins_P_Q_FILTERED.real(), "sum_sins_P_Q_FILTERED.dat");

    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(unset key)");

    // 1 step: increasing frequency in P times
    //gp.sendLine(R"(set multiplot layout 4, 1)");
    //gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_P.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_P_SPECTRA.dat" with lines, "freq_Step_1_CUTTED.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_P_FILTERED.dat" with lines)");

    // 2 step: decreasing frequency in Q times
    gp.sendLine(R"(set multiplot layout 4, 1)");
    gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_P_FILTERED.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_P_FILTERED_SPECTRA.dat" with lines, "freq_Step_2_CUTTED.dat" with lines)");
    gp.sendLine(R"(plot "sum_sins_P_Q_RESULT.dat" with lines)");

    // Results
    //gp.sendLine(R"(set multiplot layout 2, 1)");
    //gp.sendLine(R"(plot "sum_sins.dat" with lines)");
    //gp.sendLine(R"(plot "sum_sins_P_Q_RESULT.dat" with lines)");

    return 0;
}