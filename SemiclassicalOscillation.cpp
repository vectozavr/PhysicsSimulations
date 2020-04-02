//
// Created by ivan- on 31.03.2020.
//
#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

#define X_MIN 0.0001
#define X_ERROR 0.005
#define E_ERROR 0.0005

double V(double x) {
    return 4*(pow(x, -12) - pow(x, -6));
}

double x_in(double En) {
    double d_dx = 0.1;
    double d_Xin = X_MIN;

    while(d_dx > X_ERROR) {
        if( V(d_Xin) >= En ) {
            d_Xin += d_dx;
            d_dx /= 2;
        }
        d_Xin -= d_dx;
    }
    return d_Xin;
}

double x_out(double En) {
    double d_dx = 0.1;
    double d_Xout = X_MIN;

    while(d_dx > X_ERROR) {
        if( V(d_Xout) >= En ) {
            d_Xout -= d_dx;
            d_dx /= 2;
        }
        d_Xout += d_dx;
    }
    return d_Xout;
}

double testFunc(double y) {
    return y;
}

double integral_fun(double En, double x) {
    return sqrt(En - V(x));
}

double s(double gamma, double E, double En, int n, double h = 0.01) {
    double d_from = x_in(E);
    double d_to   = x_out(E);

    double d_s = 0;
    double d_x = d_from;
    while(d_x < d_to) {
        d_s += integral_fun(En, d_x) + 4*integral_fun(En, d_x + h) + integral_fun(En, d_x + 2*h);
        d_x += 2*h;
    }
    return d_s * h / 3 - (n + 0.5)*2*PI;
}

double En (int n, double E, double gamma) {
    double d_En_old = E_ERROR/2;
    double d_En = E_ERROR;
    double d_s;
    do {
        d_s = s(gamma, E, d_En, n);
        d_En -= d_s*(d_En - d_En_old)/(d_s - s(gamma, E, d_En_old, n));

    } while (d_s > E_ERROR);

}

int main() {
    double gamma = 10;
    double E = 10;

    ComplexPlot specter;

    for(int n = 0; n < 100; n++)
        specter.push(n, {En(n, E, gamma), 0});

    saveVectorPoint2DToFile(specter.real(), "specter.dat");

    GnuplotPipe gp;
    gp.sendLine(R"(plot "specter.dat" with lines)");
}

