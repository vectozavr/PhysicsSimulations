//
// Created by ivan- on 08.04.2020.
//

//
// Created by ivan- on 31.03.2020.
//
#include <vector>
#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

double runge_Kutta( double f(double, double, double, double, double),
                    double t, double p1, double p2, double p3, double p4, double h = 0.01) {
    double k1 = f(t, p1, p2, p3, p4);
    double k2 = f(t + (double)h/2, p1 + (double)k1*h/2, p2 + (double)k1*h/2, p3 + (double)k1*h/2, p4 + (double)k1*h/2);
    double k3 = f(t + (double)h/2, p1 + (double)k2*h/2, p2 + (double)k2*h/2, p3 + (double)k2*h/2, p4 + (double)k2*h/2);
    double k4 = f(t + h, p1 + h*k3, p2 + h*k3, p3 + h*k3, p4 + h*k3);

    return h*(k1 + 2*k2 + 2*k3 + k4)/6;
}

void solve_4_equation(  ComplexPlot& p1, ComplexPlot& p2, ComplexPlot& p3, ComplexPlot& p4,
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
        // Runge - Kutta method
        double diff1 = runge_Kutta(f1, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        double diff2 = runge_Kutta(f2, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        double diff3 = runge_Kutta(f3, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        double diff4 = runge_Kutta(f4, x, p1.v_c[i-1].second.real(), p2.v_c[i-1].second.real(), p3.v_c[i-1].second.real(), p4.v_c[i-1].second.real(), h);
        // ************************* //
        p1.push(x, p1.v_c[i-1].second.real() + diff1);
        p2.push(x, p2.v_c[i-1].second.real() + diff2);
        p3.push(x, p3.v_c[i-1].second.real() + diff3);
        p4.push(x, p4.v_c[i-1].second.real() + diff4);
    }
};

double testDer(double t, double p1, double p2, double p3, double p4) {
    return 1;
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

int main() {
    ComplexPlot p1;
    ComplexPlot p2;
    ComplexPlot p3;
    ComplexPlot p4;

    double x_0 = 0;
    double y_0 = 10;
    double Vx_0 = sqrt(1/y_0); //sqrt(1/y_0)
    double Vy_0 = 0;
    double t_0 = 0;
    double t_N = 200;

    double h = 0.1;

    solve_4_equation(p1, p2, p3, p4, &Vx, &Vy, &Ax, &Ay, t_0, t_N, x_0, y_0, Vx_0, Vy_0, h);

    saveVectorPoint2DToFile(p1.real(), "p1.dat");
    saveVectorPoint2DToFile(p2.real(), "p2.dat");
    saveVectorPoint2DToFile(p3.real(), "p3.dat");
    saveVectorPoint2DToFile(p4.real(), "p4.dat");

    ComplexPlot x_y;
    double error = 0;
    for(int i = 0; i < p1.size(); i += 100) {
        x_y.push(p1.v_c[i].second.real(), p2.v_c[i].second.real());
    }
    saveVectorPoint2DToFile(x_y.real(), "x_y.dat");

    GnuplotPipe gp;

    //gp.sendLine(R"(plot "p1.dat" with lines, "p2.dat" with lines, "p3.dat" with lines, "p4.dat" with lines)");
    gp.sendLine(R"(set xrange [-15:15])");
    gp.sendLine(R"(set yrange [-15:15])");
    gp.sendLine(R"(plot "x_y.dat" with lines)");

    cout << y_0 - sqrt(p1.v_c[p1.v_c.size()-1].second.real()*p1.v_c[p1.v_c.size()-1].second.real() + p2.v_c[p1.v_c.size()-1].second.real()*p2.v_c[p1.v_c.size()-1].second.real()) << endl;
}