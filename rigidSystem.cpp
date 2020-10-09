//
// Created by Иван Ильин on 08.10.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

void solve_2_equation_EX(  ComplexPlot& p1, ComplexPlot& p2,
                          double f1(double, double, double),
                          double f2(double, double, double),
                          double from, double to,
                          double p1_0, double p2_0, double h = 0.01) {
    int N = (to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's explicit method
        double diff1 = f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real())*h;
        double diff2 = f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real())*h;
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
    }
}

void solve_2_equation_IM(  ComplexPlot& p1, ComplexPlot& p2,
                           double f1(double, double, double),
                           double f2(double, double, double),
                           double from, double to,
                           double p1_0, double p2_0, double h = 0.01) {
    int N = (to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);

    for(int i = 1; i <= N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's implicit method
        double p1_2 = f1(p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), h);
        double p2_2 = f2(p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), h);
        //double p2_2 = (p1.v_c[i-1].second.real() + p2.v_c[i-1].second.real())/(1.0f + h) - p1_2;
        // ************************* //
        p1.push(x, p1_2);
        p2.push(x, p2_2);
    }
}

double du(double t, double u, double v) {
    return 998.0f*u + 1998.0f*v;
}
double dv(double t, double u, double v) {
    return -999.0f*u - 1999.0f*v;
}

double un2(double u, double v, double h) {
    return u / (1.0f + 1000.0f*h) + 1998.0f * h * (u + v) / ((1.0f + h)*(1.0f + 1000.0f*h));
}
double vn2(double u, double v, double h) {
    return (u + v) / (1 + h) - un2(u, v, h);
}

double solution(ComplexPlot& u_A, ComplexPlot& v_A, double t_0, double t_N, double alpha, double betta, double h) {
    int N = (t_N - t_0)/h;

    for(int i = 0; i < N; i++) {
        double t = t_0 + i*h;
        u_A.push(t, 2.0f*alpha*exp(-t) + betta*exp(-1000.0f*t));
        v_A.push(t, -1.0f*alpha*exp(-t) - 1.0f*betta*exp(-1000.0f*t));
    }
    saveVectorPoint2DToFile(u_A.real(), "u_A.dat");
    saveVectorPoint2DToFile(v_A.real(), "v_A.dat");
}

int main() {
    // analytical
    ComplexPlot u_A;
    ComplexPlot v_A;
    // im and ex
    ComplexPlot u_EX; ComplexPlot u_IM;
    ComplexPlot v_EX; ComplexPlot v_IM;
    // errors
    ComplexPlot u_EX_err; ComplexPlot u_IM_err;
    ComplexPlot v_EX_err; ComplexPlot v_IM_err;

    double alpha = 1;
    double betta = 5000;

    double u_0 = 2.0f*alpha + betta;
    double v_0 = -(alpha + betta);
    double t_0 = 0;
    double t_N = 0.3;

    double h = 0.00001; // Шаг интегрирования

    solution(u_A, v_A, t_0, t_N, alpha, betta, h);

    // ************************* //
    // Euler's method
    solve_2_equation_EX(u_EX, v_EX, &du, &dv, t_0, t_N, u_0, v_0, h);
    saveVectorPoint2DToFile(u_EX.real(), "u_EX.dat");
    saveVectorPoint2DToFile(v_EX.real(), "v_EX.dat");
    // Runge - Kutta method
    solve_2_equation_IM(u_IM, v_IM, &un2, &vn2, t_0, t_N, u_0, v_0, h);
    saveVectorPoint2DToFile(u_IM.real(), "u_IM.dat");
    saveVectorPoint2DToFile(v_IM.real(), "v_IM.dat");

    for(int i = 0; i < u_A.size(); i++) {
        u_EX_err.push(i*h, abs(u_A.v_c[i].second - u_EX.v_c[i].second));
        v_EX_err.push(i*h, abs(v_A.v_c[i].second - v_EX.v_c[i].second));

        u_IM_err.push(i*h, abs(u_A.v_c[i].second - u_IM.v_c[i].second));
        v_IM_err.push(i*h, abs(u_A.v_c[i].second - v_IM.v_c[i].second));
    }

    saveVectorPoint2DToFile(u_EX_err.real(), "u_EX_err.dat");
    saveVectorPoint2DToFile(v_EX_err.real(), "v_EX_err.dat");
    saveVectorPoint2DToFile(u_IM_err.real(), "u_IM_err.dat");
    saveVectorPoint2DToFile(v_IM_err.real(), "v_IM_err.dat");

    GnuplotPipe gp;

    //gp.sendLine(R"(set xrange [-15:15])");
    //gp.sendLine(R"(set yrange [-10:10])");
    //gp.sendLine(R"(plot "u_A.dat" with lines, "v_A.dat" with lines, "u_EX.dat" with lines, "v_EX.dat" with lines, "u_IM.dat" with lines, "v_IM.dat" with lines)");
    gp.sendLine(R"(plot "u_IM_err.dat" with lines, "v_IM_err.dat" with lines)");

    /*
    string plot_text = "plot ";
    for(int s = 0; s < 50; s++) {
        ComplexPlot p1_R_params;
        ComplexPlot p2_R_params;

        solve_2_equation_R2(p1_R_params, p2_R_params, &dx, &dy, t_0, t_N, u_0 + s * 0.1f, v_0 + s * 0.1f, h);
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
     */

    //gp.sendLine(R"(plot "p1_E.dat" with lines, "p2_E.dat" with lines, "p1_R2.dat" with lines, "p2_R2.dat" with lines)");
}