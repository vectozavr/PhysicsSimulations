//
// Created by Иван Ильин on 01.10.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

void solve_2_equation_R2(  ComplexPlot& p1, ComplexPlot& p2,    double f1(double, double, double),
                                                                double f2(double, double, double), double from, double to,
                                                                double p1_0, double p2_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;

        // Runge - Kutta method
        //vector<double> diff = runge_Kutta(f1, f2, f3, f4, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);

        double diff1 = h * (f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()) +
                            f1(x + h,   p1.v_c[i-1].second.real() + h*f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()),
                                        p2.v_c[i-1].second.real() + h*f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()))
                            ) / 2.0f;
        double diff2 = h * (f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()) +
                            f2(x + h,   p1.v_c[i-1].second.real() + h*f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()),
                               p2.v_c[i-1].second.real() + h*f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real()))
                            ) / 2.0f;
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
    }
}

void solve_2_equation_E(  ComplexPlot& p1, ComplexPlot& p2,
                          double f1(double, double, double),
                          double f2(double, double, double),
                          double from, double to,
                          double p1_0, double p2_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        double diff1 = f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real())*h;
        double diff2 = f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real())*h;
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
    }
}

double a_x(double t) {
    return 1.5*exp(2*t) - t*t - t - 0.5;
}

double dx(double t, double p1, double p2) {
    return 10*p1 - 2*p1*p2;
}
double dy(double t, double p1, double p2) {
    return 2*p1*p2 - 10*p2;
}

int main() {

    // MOTION OF PLANETS
    ComplexPlot p_A;

    ComplexPlot p1_R; ComplexPlot p1_E;
    ComplexPlot p2_R; ComplexPlot p2_E;

    double x_0 = 0.2;
    double y_0 = 5;
    double t_0 = 0;
    double t_N = 3;

    double h = 0.001; // Шаг интегрирования

    // ************************* //
    // Euler's method
    solve_2_equation_E(p1_E, p2_E, &dx, &dy, t_0, t_N, x_0, y_0, h);
    saveVectorPoint2DToFile(p1_E.real(), "p1_E.dat");
    saveVectorPoint2DToFile(p2_E.real(), "p2_E.dat");
    // Runge - Kutta method
    solve_2_equation_R2(p1_R, p2_R, &dx, &dy, t_0, t_N, x_0, y_0, h);
    saveVectorPoint2DToFile(p1_R.real(), "p1_R2.dat");
    saveVectorPoint2DToFile(p2_R.real(), "p2_R2.dat");
    // ************************* //

    ComplexPlot x_y_E;
    ComplexPlot x_y_R;

    for(int i = 0; i < p1_E.size(); i += 1) {
        x_y_E.push(p1_E.v_c[i].second.real(), p2_E.v_c[i].second.real());
        x_y_R.push(p1_R.v_c[i].second.real(), p2_R.v_c[i].second.real());
    }
    saveVectorPoint2DToFile(x_y_E.real(), "x_y_E.dat");
    saveVectorPoint2DToFile(x_y_R.real(), "x_y_R.dat");

    GnuplotPipe gp;

    //gp.sendLine(R"(plot "x_y_E.dat" with lines, "x_y_R.dat" with lines)");

    gp.sendLine(R"(plot "x_y_R.dat" with lines)");
    string plot_text = "plot ";
    for(int s = 0; s < 50; s++) {
        ComplexPlot p1_R_params;
        ComplexPlot p2_R_params;

        solve_2_equation_R2(p1_R_params, p2_R_params, &dx, &dy, t_0, t_N, x_0 + s * 0.1f, y_0 + s * 0.1f, h);
        saveVectorPoint2DToFile(p1_R_params.real(), to_string(s) + "_p1_R2_params.dat");
        saveVectorPoint2DToFile(p2_R_params.real(), to_string(s) + "_p2_R2_params.dat");

        ComplexPlot x_y_R_params;

        for(int j = 0; j < p1_R_params.size(); j++) {
            x_y_R_params.push(p1_R_params.v_c[j].second.real(), p2_R_params.v_c[j].second.real());
        }
        saveVectorPoint2DToFile(x_y_R_params.real(), to_string(s) + "_x_y_R_params.dat");
        plot_text += R"(")";
        plot_text += to_string(s);
        plot_text += R"(_x_y_R_params.dat" with lines notitle, )";
    }
    //gp.sendLine(plot_text);

    //gp.sendLine(R"(plot "p1_E.dat" with lines, "p2_E.dat" with lines, "p1_R2.dat" with lines, "p2_R2.dat" with lines)");
}