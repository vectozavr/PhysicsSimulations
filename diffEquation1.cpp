//
// Created by ivan- on 09.04.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

double runge_Kutta( double f(double, double), double t, double p, double h = 0.01) {
    double k1 = f(t, p);
    double k2 = f(t + (double)h/2, p + (double)k1*h/2);
    double k3 = f(t + (double)h/2, p + (double)k2*h/2);
    double k4 = f(t + h, p + h*k3);

    return h*(k1 + 2*k2 + 2*k3 + k4)/6;
}

void solve_1_equation_R(  ComplexPlot& p, double f(double, double), double from, double to, double p_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p.push(from, p_0);

    for(int i = 1; i < N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        //double diff = f(x, p.v_c[i-1].second.real())*h;
        // ************************* //
        // Runge - Kutta method
        double diff = runge_Kutta(f, x, p.v_c[i-1].second.real(), h);
        // ************************* //
        p.push(x, p.v_c[i-1].second.real() + diff);
    }
}

void solve_1_equation_E(  ComplexPlot& p, double f(double, double), double from, double to, double p_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p.push(from, p_0);

    for(int i = 1; i < N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        double diff = f(x, p.v_c[i-1].second.real())*h;
        // ************************* //
        // Runge - Kutta method
        //double diff = runge_Kutta(f, x, p.v_c[i-1].second.real(), h);
        // ************************* //
        p.push(x, p.v_c[i-1].second.real() + diff);
    }
}

double fun(double t, double p) {
    return 2*(t*t + p);
}

double analytical(double t) {
    return 1.5*exp(2*t) - t*t - t - 0.5;
}

int main() {
    ComplexPlot p_R;
    ComplexPlot p_E;

    double y_0 = 1;

    double t_0 = 0;
    double t_N = 1;

    double h = 0.1;

    // ************************* //
    // Euler's method
    solve_1_equation_E(p_E, &fun, t_0, t_N, y_0, h);
    saveVectorPoint2DToFile(p_E.real(), "p_E.dat");
    // Runge - Kutta method
    solve_1_equation_R(p_R, &fun, t_0, t_N, y_0, h);
    saveVectorPoint2DToFile(p_R.real(), "p_R.dat");
    // ************************* //

    ComplexPlot p_analytical;
    for(int i = 0; i < p_R.size(); i++)
        p_analytical.push(p_R.v_c[i].first, analytical(p_R.v_c[i].first));
    saveVectorPoint2DToFile(p_R.real(), "p_analytical.dat");

    GnuplotPipe gp;

    double d_analytical = p_analytical  .v_c.back().second.real();
    double d_E          = p_E           .v_c.back().second.real();
    double d_R          = p_R           .v_c.back().second.real();
    cout << "analytical y = "   << d_analytical << endl;
    cout << "Euler y = "        << d_E          << endl;
    cout << "Runge y = "        << d_R          << endl;
    cout << "h = "              << h            << endl;

    //gp.sendLine(R"(set xrange [-15:15])");
    //gp.sendLine(R"(set yrange [-15:15])");set label 2 "A = %g", a at 12,1200
    gp.sendLine(R"(plot "p_R.dat" with lines, "p_E.dat" with lines, "p_analytical.dat" with lines)");
}