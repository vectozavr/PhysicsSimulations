//
// Created by Иван Ильин on 25.09.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

double bessel(double x, int m, double h) {
    double t = 0.0f;
    double S = 0.0f;

    while (t < PI) {
        S += cos(m*t - x*sin(t)) + 4.0f*cos(m*(t + h/2.0f) - x*sin(t + h/2.0f)) + cos(m*(t + h) - x*sin(t + h));
        t += h;
    }

    return S * h / 6.0f / PI;
}

double bessel0(double x, double h) {
    return bessel(x, 0, h);
}

double bessel1(double x, double h) {
    return bessel(x, 1, h);
}

double dbessel0(double x, double h) {
    return (bessel0(x+h, h) - bessel0(x-h, h))/(2.0f*h);
}

int main() {

    ComplexPlot errors_with_diff_h;

    double h = 0.000001f;
    int N = 10;


    while (h <= 1) {

        ComplexPlot errors;
        for(int i = 0; i < N; i++) {
            double x = 2.0f * PI * i / N;
            errors.push(x, {abs(dbessel0(x, h) + bessel1(x, h)),0});
        }

        errors_with_diff_h.push(log(h), {log(errors.sumReal()), 0});

        h *= 10;
    }

    saveVectorPoint2DToFile(errors_with_diff_h.real(), "errors_with_diff_h.dat");

    // GRAPH PLOT
    GnuplotPipe gp;
    //gp.sendLine(R"(set multiplot layout 2, 2)");

    gp.sendLine(R"(plot "errors_with_diff_h.dat" with lines)");
}
