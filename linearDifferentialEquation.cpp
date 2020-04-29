//
// Created by ivan- on 08.04.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

vector<double> runge_Kutta( double f1(double, double, double, double, double),
                    double f2(double, double, double, double, double),
                    double f3(double, double, double, double, double),
                    double f4(double, double, double, double, double),
                    double t, double p1, double p2, double p3, double p4, double h = 0.01) {
    vector<double> result;

    double k1_1 = f1(t, p1, p2, p3, p4);
    double k1_2 = f2(t, p1, p2, p3, p4);
    double k1_3 = f3(t, p1, p2, p3, p4);
    double k1_4 = f4(t, p1, p2, p3, p4);

    double k2_1 = f1(t + h/2, p1 + k1_1*h/2, p2 + k1_2*h/2, p3 + k1_3*h/2, p4 + k1_4*h/2);
    double k2_2 = f2(t + h/2, p1 + k1_1*h/2, p2 + k1_2*h/2, p3 + k1_3*h/2, p4 + k1_4*h/2);
    double k2_3 = f3(t + h/2, p1 + k1_1*h/2, p2 + k1_2*h/2, p3 + k1_3*h/2, p4 + k1_4*h/2);
    double k2_4 = f4(t + h/2, p1 + k1_1*h/2, p2 + k1_2*h/2, p3 + k1_3*h/2, p4 + k1_4*h/2);

    double k3_1 = f1(t + h/2, p1 + k2_1*h/2, p2 + k2_2*h/2, p3 + k2_3*h/2, p4 + k2_4*h/2);
    double k3_2 = f2(t + h/2, p1 + k2_1*h/2, p2 + k2_2*h/2, p3 + k2_3*h/2, p4 + k2_4*h/2);
    double k3_3 = f3(t + h/2, p1 + k2_1*h/2, p2 + k2_2*h/2, p3 + k2_3*h/2, p4 + k2_4*h/2);
    double k3_4 = f4(t + h/2, p1 + k2_1*h/2, p2 + k2_2*h/2, p3 + k2_3*h/2, p4 + k2_4*h/2);

    double k4_1 = f1(t + h, p1 + k3_1*h, p2 + k3_2*h, p3 + k3_3*h, p4 + k3_4*h);
    double k4_2 = f2(t + h, p1 + k3_1*h, p2 + k3_2*h, p3 + k3_3*h, p4 + k3_4*h);
    double k4_3 = f3(t + h, p1 + k3_1*h, p2 + k3_2*h, p3 + k3_3*h, p4 + k3_4*h);
    double k4_4 = f4(t + h, p1 + k3_1*h, p2 + k3_2*h, p3 + k3_3*h, p4 + k3_4*h);

    result.push_back(h*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6);
    result.push_back(h*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6);
    result.push_back(h*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3)/6);
    result.push_back(h*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4)/6);

    return result;
}

void solve_4_equation_R(  ComplexPlot& p1, ComplexPlot& p2, ComplexPlot& p3, ComplexPlot& p4,
                        double f1(double, double, double, double, double),
                        double f2(double, double, double, double, double),
                        double f3(double, double, double, double, double),
                        double f4(double, double, double, double, double),
                        double from, double to,
                        double p1_0, double p2_0, double p3_0, double p4_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);
    p3.push(from, p3_0);
    p4.push(from, p4_0);

    for(int i = 1; i < N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        //double diff1 = f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        //double diff2 = f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        //double diff3 = f3(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        //double diff4 = f4(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        // ************************* //
        // Runge - Kutta method (OLD)
        vector<double> diff = runge_Kutta(f1, f2, f3, f4, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);

        double diff1 = diff[0];
        double diff2 = diff[1];
        double diff3 = diff[2];
        double diff4 = diff[3];
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
        p3.push(x, p3.v_c[i-1].second.real() + diff3);
        p4.push(x, p4.v_c[i-1].second.real() + diff4);
    }
}

void solve_4_equation_E(  ComplexPlot& p1, ComplexPlot& p2, ComplexPlot& p3, ComplexPlot& p4,
                          double f1(double, double, double, double, double),
                          double f2(double, double, double, double, double),
                          double f3(double, double, double, double, double),
                          double f4(double, double, double, double, double),
                          double from, double to,
                          double p1_0, double p2_0, double p3_0, double p4_0, double h = 0.01) {
    int N = (int)(to - from)/h;
    p1.push(from, p1_0);
    p2.push(from, p2_0);
    p3.push(from, p3_0);
    p4.push(from, p4_0);

    for(int i = 1; i < N; i++) {
        double x = from + i*h;
        // ************************* //
        // Euler's method
        double diff1 = f1(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        double diff2 = f2(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        double diff3 = f3(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        double diff4 = f4(x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real())*h;
        // ************************* //
        // Runge - Kutta method
        //double diff1 = runge_Kutta(f1, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        //double diff2 = runge_Kutta(f2, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        //double diff3 = runge_Kutta(f3, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        //double diff4 = runge_Kutta(f4, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
        p3.push(x, p3.v_c[i-1].second.real() + diff3);
        p4.push(x, p4.v_c[i-1].second.real() + diff4);
    }
}

double fun(double t, double p1, double p2, double p3, double p4) {
    return 2*(t*t + p1);
}

double analytical(double t) {
    return 1.5*exp(2*t) - t*t - t - 0.5;
}

double Vx(double t, double p1, double p2, double p3, double p4) {
    return p3;
}
double Vy(double t, double p1, double p2, double p3, double p4) {
    return p4;
}
double Ax(double t, double p1, double p2, double p3, double p4) {
    return -p1/pow((p1*p1 + p2*p2), 1.5);
}
double Ay(double t, double p1, double p2, double p3, double p4) {
    return -p2/pow((p1*p1 + p2*p2), 1.5);
}

void circle(ComplexPlot& p, double R, int N = 100) {
    for(int i = 0; i <= N; i++) {
        double alpha = i*2*PI/N;
        double x = R*cos(alpha);
        double y = R*sin(alpha);
        p.push(x, y);
    }
}

int main() {
    // TEST OF FIRST PROGRAM
    /*
    ComplexPlot p_NULL;

    ComplexPlot p_R;
    ComplexPlot p_E;

    double y_0 = 1;

    double t_0 = 0;
    double t_N = 1;

    double h = 0.1;

    // ************************* //
    // Euler's method
    solve_4_equation_E(p_E, p_NULL, p_NULL, p_NULL, &fun, &fun, &fun, &fun, t_0, t_N, y_0, 0, 0, 0, h);
    saveVectorPoint2DToFile(p_E.real(), "p_E.dat");
    // Runge - Kutta method
    solve_4_equation_R(p_R, p_NULL, p_NULL, p_NULL, &fun, &fun, &fun, &fun, t_0, t_N, y_0, 0, 0, 0, h);
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

    gp.sendLine(R"(plot "p_R.dat" with lines, "p_E.dat" with lines, "p_analytical.dat" with lines)");
    */
    // MOTION OF PLANETS
    ComplexPlot circl;

    ComplexPlot p1_R; ComplexPlot p1_E;
    ComplexPlot p2_R; ComplexPlot p2_E;
    ComplexPlot p3_R; ComplexPlot p3_E;
    ComplexPlot p4_R; ComplexPlot p4_E;

    double x_0 = 0;     // Начальная позиция планеты по X
    double y_0 = 10;    // Начальная позиция планеты по Y
    double Vx_0 = sqrt(1/y_0);  // <- Скорость Vx_0 такая,
    double Vy_0 = 0;                // что поддерживается
    double t_0 = 0; // Время начала // центростремительное
    double t_N = 200;// Время окончания // ускорение.

    double h = 10; // Шаг интегрирования

    // ************************* //
    // Euler's method
    solve_4_equation_E(p1_E, p2_E, p3_E, p4_E, &Vx, &Vy, &Ax, &Ay, t_0, t_N, x_0, y_0, Vx_0, Vy_0, h);
    saveVectorPoint2DToFile(p1_E.real(), "p1_E.dat");
    saveVectorPoint2DToFile(p2_E.real(), "p2_E.dat");
    saveVectorPoint2DToFile(p3_E.real(), "p3_E.dat");
    saveVectorPoint2DToFile(p4_E.real(), "p4_E.dat");
    // Runge - Kutta method
    solve_4_equation_R(p1_R, p2_R, p3_R, p4_R, &Vx, &Vy, &Ax, &Ay, t_0, t_N, x_0, y_0, Vx_0, Vy_0, h);
    saveVectorPoint2DToFile(p1_R.real(), "p1_R.dat");
    saveVectorPoint2DToFile(p2_R.real(), "p2_R.dat");
    saveVectorPoint2DToFile(p3_R.real(), "p3_R.dat");
    saveVectorPoint2DToFile(p4_R.real(), "p4_R.dat");
    // ************************* //
    circle(circl, y_0);
    saveVectorPoint2DToFile(circl.real(), "circle.dat");

    ComplexPlot x_y_E;
    ComplexPlot x_y_R;
    double error = 0;
    for(int i = 0; i < p1_R.size(); i += 1) {
        x_y_E.push(p1_E.v_c[i].second.real(), p2_E.v_c[i].second.real());
        x_y_R.push(p1_R.v_c[i].second.real(), p2_R.v_c[i].second.real());
    }
    saveVectorPoint2DToFile(x_y_E.real(), "x_y_E.dat");
    saveVectorPoint2DToFile(x_y_R.real(), "x_y_R.dat");

    GnuplotPipe gp;

    double d_analytical = sqrt(circl.v_c.back().second.real()*circl.v_c.back().second.real() + circl.v_c.back().first*circl.v_c.back().first);
    double d_E          = sqrt(x_y_E.v_c.back().second.real()*x_y_E.v_c.back().second.real() + x_y_E.v_c.back().first*x_y_E.v_c.back().first);
    double d_R          = sqrt(x_y_R.v_c.back().second.real()*x_y_R.v_c.back().second.real() + x_y_R.v_c.back().first*x_y_R.v_c.back().first);
    cout << "analytical R = "   << d_analytical << endl;
    cout << "Euler R = "        << d_E          << endl;
    cout << "Runge R = "        << d_R          << endl;
    cout << "h = "              << h            << endl;

    //gp.sendLine(R"(plot "p1.dat" with lines, "p2.dat" with lines, "p3.dat" with lines, "p4.dat" with lines)");
    gp.sendLine(R"(set xrange [-15:15])");
    gp.sendLine(R"(set yrange [-15:15])");
    gp.sendLine(R"(plot "x_y_E.dat" with lines, "x_y_R.dat" with lines, "circle.dat" with lines)");
}