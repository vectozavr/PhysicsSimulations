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

#define X_MIN pow(2, (double)1/6)
#define X_ERROR 0.0005
#define E_ERROR 0.0005

double V(double x) {
    return 4*(pow(x, -12) - pow(x, -6));
}

double x_in(double En) {
    double d_dx = 0.1;
    double d_Xin = X_MIN;

    while(d_dx > X_ERROR && (d_Xin < 10)) {
        if( V(d_Xin) >= En ) {
            d_Xin += d_dx;
        } else {
            d_Xin -= d_dx;
            d_dx /= 2;
        }
    }
    return d_Xin;
}

double x_out(double En) {
    double d_dx = 0.1;
    double d_Xout = 3;

    while(d_dx > X_ERROR) {
        if( V(d_Xout) >= En && (d_Xout > 0)) {
            d_Xout -= d_dx;
        } else {
            d_Xout += d_dx;
            d_dx /= 2;
        }
    }
    return d_Xout;
}

double integral_fun(double En, double x) {

    if(En - V(x) > 0)
        return sqrt(En - V(x));
    else
        return 0;
}

double s(double E, double gamma, double h = 0.01) {
    double d_from = x_in(E);
    double d_to   = x_out(E);

    double d_s = 0;
    double d_x = d_from;
    while(d_x < d_to) {
        d_s += integral_fun(E, d_x) + 4*integral_fun(E, d_x + h) + integral_fun(E, d_x + 2*h);
        d_x += 2*h;
    }
    return gamma * d_s * h / 3;
}

int N(double E, double gamma, double h = 0.01) {
    return (int)(s(E, gamma, h)/PI - 0.5);
}

double En (double gamma, int n, double h = 0.01) {
    double E1 = -1;
    double E2 = E1 + abs(E1)/4;
    double F1 = -(double)PI/2;
    double F2 = 0;

    double dE = 2*E_ERROR;
    double E = 0;

    while (abs(dE) > E_ERROR) {
        E = E2;
        F2 = s(E, gamma, h) - (n + 0.5)*PI;
        if(F2 == F1)
            break;
        dE = -F2*(E2 - E1)/(F2 - F1);
        E1 = E2;
        F1 = F2;
        E2 = E1 + dE;
        if(E2 > 0) E2 = -E_ERROR;
    }

    return E;
}

int main() {
    double gamma = 100;

    ComplexPlot potentialV;
    for(int k = 10; k < 1000; k++) {
        potentialV.push((double)k/100, {V((double)k/100), 0});
    }
    saveVectorPoint2DToFile(potentialV.real(), "p.dat");

    vector<ComplexPlot> plots;

    int i_N = N(-E_ERROR, gamma); // Кол-во уровней.
    for(int i = 0; i < i_N; i++) {
        double E = En(gamma, i);
        ComplexPlot plot;

        plot.push(x_in(E), {E, 0});
        plot.push(x_out(E), {E, 0});

        string name = "plot_" + to_string(i);

        saveVectorPoint2DToFile(plot.real(), name);
    }

    GnuplotPipe gp;
    gp.sendLine(R"(set yrange [-1:1])");
    gp.sendLine(R"(set xrange [0:3])");

    string gpLine;
    gpLine += R"(plot "p.dat" with lines, )";
    for(int k = 0; k < i_N; k++) {
        string name = "plot_" + to_string(k);
        gpLine += R"(")";
        gpLine += name;
        gpLine += R"("with lines,)";
    }
    gp.sendLine(gpLine);

    cout << N(-E_ERROR, gamma) << endl;
    //gp.sendLine(R"(plot "specter1.dat" with lines,  "specter2.dat" with lines)");
}

