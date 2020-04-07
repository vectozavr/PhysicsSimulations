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

#define X_MIN 0.5
#define X_ERROR 0.005
#define E_ERROR 0.9

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

double s(double gamma, double E, double En, int n, double h = 0.01) {
    double d_from = x_in(E);
    double d_to   = x_out(E);

    double d_s = 0;
    double d_x = d_from;
    while(d_x < d_to) {
        d_s += integral_fun(En, d_x) + 4*integral_fun(En, d_x + h) + integral_fun(En, d_x + 2*h);
        d_x += 2*h;
    }
    return gamma * d_s * h / 3 - (n + 0.5)*2*PI;
}

double En (int n, double E, double gamma) {
    double d_En_old = -1;
    double d_En = d_En_old + abs(d_En_old)/4;
    double d_s;
    do {
        d_s = s(gamma, E, d_En, n);
        d_En -= d_s*(d_En - d_En_old)/(d_s - s(gamma, E, d_En_old, n));
    } while (abs(d_s) > E_ERROR);
}

int main() {
    double gamma = 10;
    double E1 = 10;
    double E2 = 20;


    ComplexPlot potentialV;
    for(int k = 10; k < 1000; k++) {
        potentialV.push((double)k/100, {V((double)k/100), 0});
    }
    saveVectorPoint2DToFile(potentialV.real(), "p.dat");

    vector<ComplexPlot> plots;
    int d_plots = 5;
    for(int i = 0; i < d_plots; i++) {
        double E = En(i*3+10, -0.1, gamma);
        //double E = -(double)(i+1)/10;
        ComplexPlot plot;

        plot.push(x_in(-E), {-E, 0});
        plot.push(x_out(-E), {-E, 0});

        string name = "plot_" + to_string(i);

        saveVectorPoint2DToFile(plot.real(), name);
    }

    GnuplotPipe gp;
    gp.sendLine(R"(set yrange [-1:1])");
    gp.sendLine(R"(set xrange [0:3])");

    string gpLine;
    gpLine += R"(plot "p.dat" with lines, )";
    for(int k = 0; k < d_plots; k++) {
        string name = "plot_" + to_string(k);
        gpLine += R"(")";
        gpLine += name;
        gpLine += R"("with lines,)";
    }
    gp.sendLine(gpLine);

    cout << x_in(-0.2) << endl;
    cout << x_in(-0.7) << endl;
    //gp.sendLine(R"(plot "specter1.dat" with lines,  "specter2.dat" with lines)");
}

