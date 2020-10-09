//
// Created by Иван Ильин on 18.09.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

double fun1(double x) {
    return 1.0f / (1.0f + x*x);
}

double fun2(double x) {
    return pow(x, 1.0f / 3.0f) * exp(sin(x));
}

double traps(double f(double x), double x0, double x1, int N) {
    double h = (x1 - x0) / N;
    double x = x0;

    double S = 0.0f;

    while (x < x1) {
        S += f(x) + f(x+h);
        x += h;
    }

    return S * h / 2;
}

double simps(double f(double x), double x0, double x1, int N) {
    double h = (x1 - x0) / N;
    double x = x0;

    double S = 0.0f;

    while (x < x1) {
        S += f(x) + 4.0f*f(x + h/2.0f) + f(x + h);
        x += h;
    }

    return S * h / 6;
}

int main() {

    int N = 4;

    cout << "----------------------" << endl;
    cout << "Traps method" << endl;
    cout << traps(fun1, -1.0f, 1.0f, N) << endl;
    cout << traps(fun2, 0.1f, 1.0f, N) << endl;
    cout << "----------------------" << endl;
    cout << "Simps method" << endl;
    cout << simps(fun1, -1.0f, 1.0f, N) << endl;
    cout << simps(fun2, 0.1f, 1.0f, N) << endl;

    ComplexPlot error_traps_fun1;
    ComplexPlot error_traps_fun2;

    ComplexPlot error_simps_fun1;
    ComplexPlot error_simps_fun2;

    vector<pair<pair<double, double>, pair<double, double>>> errors;

    for(; N <= 128; N *= 2) {

        double trap1 = traps(fun1, -1.0f, 1.0f, N);
        double trap2 = traps(fun2, 0.1f, 1.0f, N);
        error_traps_fun1.push(N, trap1);
        error_traps_fun2.push(N, trap2);

        double simps1 = simps(fun1, -1.0f, 1.0f, N);
        double simps2 = simps(fun2, 0.1f, 1.0f, N);
        error_simps_fun1.push(N, simps1);
        error_simps_fun2.push(N, simps2);

        double first_int = PI / 2.0f;
        double second_int = 1.259003329982857;

        //cout << "-----------------" << endl;
        //cout << "N = " << N << endl;
        //cout << "Error traps f1 = " << abs(trap1 - first_int) << endl;
        //cout << "Error traps f2 = " << abs(trap2 - second_int) << endl;
        //cout << "-" << endl;
        //cout << "Error simps f1 = " << abs(simps1 - first_int) << endl;
        //cout << "Error simps f2 = " << abs(simps2 - second_int) << endl;

        errors.push_back({{abs(trap1 - first_int), abs(trap2 - second_int)}, {abs(simps1 - first_int), abs(simps2 - second_int)}});
    }

    saveVectorPoint2DToFile(error_traps_fun1.real(), "error_traps_fun1.dat");
    saveVectorPoint2DToFile(error_traps_fun2.real(), "error_traps_fun2.dat");
    saveVectorPoint2DToFile(error_simps_fun1.real(), "error_simps_fun1.dat");
    saveVectorPoint2DToFile(error_simps_fun2.real(), "error_simps_fun2.dat");


    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 2, 2)");

    gp.sendLine(R"(plot "error_traps_fun1.dat" with lines)");
    gp.sendLine(R"(plot "error_traps_fun2.dat" with lines)");
    gp.sendLine(R"(plot "error_simps_fun1.dat" with lines)");
    gp.sendLine(R"(plot "error_simps_fun2.dat" with lines)");

    for(int i = 0; i < errors.size()-1; i++) {
        //cout << "Et1(N)/Et1(2N) = " << errors[i].first.first / errors[i+1].first.first << endl;
        //cout << "Et2(N)/Et2(2N) = " << errors[i].first.second / errors[i+1].first.second << endl;
        //cout << "-" << endl;
        //cout << "Es1(N)/Es1(2N) = " << errors[i].second.first / errors[i+1].second.first << endl;
        cout << "Es2(N)/Es2(2N) = " << errors[i].second.second / errors[i+1].second.second << endl;

        cout << "-----------------" << endl;
        cout << "-----------------" << endl;
    }
}