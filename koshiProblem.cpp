//
// Created by Иван Ильин on 25.09.2020.
//
#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

double runge_Kutta4( double f(double, double), double t, double p, double h = 0.01f) {

    double k1_1 = f(t, p);
    double k2_1 = f(t + h / 2.0f, p + k1_1*h / 2.0f);
    double k3_1 = f(t + h / 2.0f, p + k2_1*h / 2.0f);
    double k4_1 = f(t + h, p + k3_1*h);

    return h * (k1_1 + 2.0f * k2_1 + 2.0f * k3_1 + k4_1) / 6.0f;
}

void solve_equation_R4(  ComplexPlot& p, double f(double, double), double from, double to, double p_0, double h = 0.01f) {
    int N = (to - from) / h;
    p.push(from, p_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;

        // Runge - Kutta 4 method
        double diff = runge_Kutta4(f, x, p.v_c[i-1].second.real(), h);

        // ************************* //
        p.push(x, p.v_c[i-1].second.real() + diff);
    }
}

void solve_equation_R2(  ComplexPlot& p, double f(double, double), double from, double to, double p_0, double h = 0.01f) {
    int N = (to - from) / h;
    p.push(from, p_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;

        // Runge - Kutta 2 method
        double diff = h * (f(x, p.v_c[i-1].second.real()) + f(x + h, p.v_c[i-1].second.real() + h*f(x, p.v_c[i-1].second.real()))) / 2.0f;

        // ************************* //
        p.push(x, p.v_c[i-1].second.real() + diff);
    }
}

void solve_equation_E(  ComplexPlot& p, double f(double, double), double from, double to, double p_0, double h = 0.01) {
    int N = (to - from)/h;
    p.push(from, p_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        double diff = f(x, p.v_c[i-1].second.real()) * h;

        // ************************* //
        p.push(x, p.v_c[i-1].second.real() + diff);
    }
}

double dx1(double t, double x) {
    return -x;
}

double x1(double t) {
    return exp(-t);
}

int main() {
    ComplexPlot p_A;

    ComplexPlot p_E;
    ComplexPlot p_R2;
    ComplexPlot p_R4;

    ComplexPlot p_E_error;
    ComplexPlot p_R2_error;
    ComplexPlot p_R4_error;

    double x_0 = 1; // x(0)
    double t_0 = 0; // Время начала
    double t_N = 3; // Время окончания

    double h = 0.01f; // Шаг интегрирования

    // ************************* //
    // Analytical solution
    for(int i = 0; i <= (t_N-t_0) / h; i++) {
        p_A.push(t_0 + i*h, x1(t_0 + i*h));
    }
    saveVectorPoint2DToFile(p_A.real(), "p_A.dat");

    // ************************* //
    // Euler's method
    solve_equation_E(p_E, &dx1, t_0, t_N, x_0, h);
    saveVectorPoint2DToFile(p_E.real(), "p_E.dat");

    // ************************* //
    // Runge Kutta 2 method
    solve_equation_R2(p_R2, &dx1, t_0, t_N, x_0, h);
    saveVectorPoint2DToFile(p_R2.real(), "p_R2.dat");
    // ************************* //
    // Runge Kutta 4 method
    solve_equation_R4(p_R4, &dx1, t_0, t_N, x_0, h);
    saveVectorPoint2DToFile(p_R4.real(), "p_R4.dat");

    // ************************* //
    // Errors
    for(int i = 0; i < p_A.size(); i++) {
        p_E_error.push(t_0 + i*h,   abs(p_E.v_c[i].second.real()     - p_A.v_c[i].second.real()));
        p_R2_error.push(t_0 + i*h,  abs(p_R2.v_c[i].second.real()    - p_A.v_c[i].second.real()));
        p_R4_error.push(t_0 + i*h,  abs(p_R4.v_c[i].second.real()    - p_A.v_c[i].second.real()));
    }
    saveVectorPoint2DToFile(p_E_error.real(), "p_E_error.dat");
    saveVectorPoint2DToFile(p_R2_error.real(), "p_R2_error.dat");
    saveVectorPoint2DToFile(p_R4_error.real(), "p_R4_error.dat");

    ComplexPlot p_errors;
    for(int i = 1; i < 6; i++) {
        double _h = pow(10, -i);
        ComplexPlot sol;
        ComplexPlot err;
        solve_equation_R4(sol, &dx1, t_0, t_N, x_0, _h);
        for(int s = 0; s < sol.size(); s++) {
            err.push(t_0 + i*h,   abs(sol.v_c[i].second.real() - p_A.v_c[i].second.real()));
        }
        p_errors.push(log(_h),log(err.midReal()));
    }
    saveVectorPoint2DToFile(p_errors.real(), "p_errors.dat");



    // GRAPH PLOT
    GnuplotPipe gp;
    //gp.sendLine(R"(set multiplot layout 3, 1)");
    //gp.sendLine(R"(plot "p_E_error.dat" with lines)");
    //gp.sendLine(R"(plot "p_R2_error.dat" with lines)");
    //gp.sendLine(R"(plot "p_R4_error.dat" with lines)");

    gp.sendLine(R"(plot "p_errors.dat" with lines)");
}