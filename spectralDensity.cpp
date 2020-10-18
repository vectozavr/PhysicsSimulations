//
// Created by Иван Ильин on 15.10.2020.
//

#include <algorithm>
#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;
using namespace std;

int main() {

    double t_0 = 0.0f;
    double t_1 = 2.0f*PI;

    double omega_0 = 5.1f;
    double omega_1 = 5.0f*omega_0;

    double a_0 = 1.0f;
    double a_1 = 0.002f;

    double tau = 0.01f;

    int N = 1.0f + (t_1 - t_0) / tau;

    ComplexPlot f;
    ComplexPlot hann_window;
    ComplexPlot cross_f_hann;

    for(int i = 0; i < N; i++) {
        double y = (a_0 * sin(omega_0 * (t_0 + i * tau)) + a_1 * sin(omega_1 * (t_0 + i * tau)));

        f.push(t_0 + i * tau, y);
        hann_window.push(t_0 + i * tau, 0.5f*(1 - cos(2.0f*PI*i/(N-1))));
    }

    cross(f, hann_window, cross_f_hann);

    //cross(sum_sins_P_FILTERED_SPECTRA, freq_Step_2_CUTTED, cross_sum_sins2);
    //inverseFourierTransform(cross_sum_sins2, sum_sins_P_Q_FILTERED);
    //fourierTransform(sum_sins_P_Q_FILTERED, sum_sins_P_Q_FILTERED_SPECTRA);

    ComplexPlot f_SPECTRA;
    ComplexPlot hann_window_SPECTRA;
    ComplexPlot cross_f_hann_SPECTRA;

    ComplexPlot f_cross;
    ComplexPlot f_cross_SPECTRA;
    crossCorrelation(f, f, f_cross);
    fourierTransform(f_cross, f_cross_SPECTRA);
    saveVectorPoint2DToFile(f_cross.real(), "f_cross.dat");
    saveVectorPoint2DToFile(f_cross_SPECTRA.abs(), "f_cross_SPECTRA.dat", f_cross_SPECTRA.size()/10);

    fourierTransform(f, f_SPECTRA);
    fourierTransform(hann_window, hann_window_SPECTRA);
    fourierTransform(cross_f_hann, cross_f_hann_SPECTRA);

    saveVectorPoint2DToFile(f.real(), "f.dat");
    saveVectorPoint2DToFile(hann_window.real(), "hann_window.dat");
    saveVectorPoint2DToFile(cross_f_hann.real(), "cross_f_hann.dat");

    for(int i = 0; i < f_SPECTRA.size(); i++) {
        double abs = sqrt(f_SPECTRA.v_c[i].second.real() * f_SPECTRA.v_c[i].second.real() +
                                  f_SPECTRA.v_c[i].second.imag() * f_SPECTRA.v_c[i].second.imag());
        f_SPECTRA.v_c[i].second = log(abs);
    }
    for(int i = 0; i < cross_f_hann_SPECTRA.size(); i++) {
        double abs = sqrt(cross_f_hann_SPECTRA.v_c[i].second.real() * cross_f_hann_SPECTRA.v_c[i].second.real() +
                                  cross_f_hann_SPECTRA.v_c[i].second.imag() * cross_f_hann_SPECTRA.v_c[i].second.imag());
        cross_f_hann_SPECTRA.v_c[i].second = log(abs);
    }

    saveVectorPoint2DToFile(f_SPECTRA.real(), "f_SPECTRA.dat", f_SPECTRA.size()/10);
    saveVectorPoint2DToFile(hann_window_SPECTRA.abs(), "hann_window_SPECTRA.dat", f_SPECTRA.size()/10);
    saveVectorPoint2DToFile(cross_f_hann_SPECTRA.real(), "cross_f_hann_SPECTRA.dat", f_SPECTRA.size()/10);

    // GRAPH PLOT
    GnuplotPipe gp;

    gp.sendLine(R"(set multiplot layout 3, 2)");

    gp.sendLine(R"(plot "f.dat" with lines)");
    gp.sendLine(R"(plot "f_SPECTRA.dat" with lines)");
    //gp.sendLine(R"(plot "f_cross.dat" with lines)");
    //gp.sendLine(R"(plot "f_cross_SPECTRA.dat" with lines)");

    gp.sendLine(R"(plot "hann_window.dat" with lines)");
    gp.sendLine(R"(plot "hann_window_SPECTRA.dat" with lines)");

    gp.sendLine(R"(plot "cross_f_hann.dat" with lines)");
    gp.sendLine(R"(plot "cross_f_hann_SPECTRA.dat" with lines)");

    return 0;
}