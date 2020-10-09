//
// Created by Иван Ильин on 25.09.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

double xk(unsigned k, unsigned n) {
    return 1.0f + (double)k / (double)n;
}

double li(double x, unsigned i, unsigned n) {
    double res = 1.0f;

    for(int k = 0; k <= n; k++)
        if(i != k)
            res *= (x - xk(k, n));

    return res;
}

double yk(unsigned k, unsigned n) {
    return log(xk(k, n));
}

double pn(double x, unsigned n) {
    double res = 0.0f;

    for(int k = 0; k <= n; k++)
        res += yk(k, n) * li(x, k, n) / li(xk(k, n), k, n);

    return res;
}

int main() {

    double x0 = 1.0f;
    double x1 = 4.0f;
    double h = 0.01f;

    //drawing log(x)
    ComplexPlot plot_log;
    double t = x0;
    while (t <= x1) {
        plot_log.push(t, {log(t), 0});
        t += h;
    }
    saveVectorPoint2DToFile(plot_log.real(), "plot_log.dat");

    ComplexPlot plot_errors;

    for(int n = 4; n <= 30; n++) {
        ComplexPlot plot;

        double x = x0;
        while (x <= x1) {
            plot.push(x, {pn(x, n) - log(x), 0});
            x += h;
        }

        plot_errors.push(n, {log(plot.maxReal()), 0});

        saveVectorPoint2DToFile(plot.real(), "plot" + to_string(n) + ".dat");
    }

    saveVectorPoint2DToFile(plot_errors.real(), "plot_errors.dat");

    // GRAPH PLOT
    GnuplotPipe gp;
    gp.sendLine(R"(plot "plot_log.dat" with lines , "plot4.dat" , "plot5.dat" , "plot6.dat" , "plot7.dat" , "plot8.dat" , "plot9.dat" , "plot10.dat" , "plot11.dat" , "plot12.dat" , "plot13.dat" , "plot14.dat" , "plot15.dat" )");

    //gp.sendLine(R"(plot "plot5.dat" with lines, "plot6.dat" with lines, "plot22.dat" with lines)");

    //gp.sendLine(R"(plot "plot_errors.dat")");

    return 0;
}