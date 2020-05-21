//
// Created by ivan- on 21.05.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"
#include "Matrix.hpp"

using namespace std;
using namespace vemath;

#define PI 3.1415926535

unsigned fact(unsigned i) {
    unsigned result = 1;
    for(int k = 2; k <= i; k++)
        result *= k;
    return result;
}

double alpha_i(double t, unsigned i) {
    return pow(-1, i) * pow(t, 2*i);
}

double betta_i(double t, unsigned i) {
    return pow(t, 2*i) / fact(i);
}

double funct(double t = 0, unsigned i = 0) {
    return 1;
}

double integral(double f1(double t, unsigned i), double f2(double t, unsigned i), double x_from, double x_to, unsigned i, unsigned j, double h = 0.01) {

    double s = 0;
    double x = x_from;
    while(x < x_to) {
        s += f1(x, i)*f2(x, j) + 4*f1(x + h, i)*f2(x + h, j) + f1(x + 2*h, i)*f2(x + 2*h, j);
        x += 2*h;
    }
    return s * h / 3;
}

int sign(double x) {
    return x > 0 ? 1 : -1;
}

void Degenerate_Fred_2(ComplexPlot& y, double x_from, double x_to, double lambda, double alpha(double t, unsigned i), double betta(double t, unsigned i), double f(double t, unsigned i),
                       unsigned n = 2, double h = 0.01) {
    Matrix r(n+1, 1);
    Matrix M(n+1, n+1);

    for(int i = 0; i <= n; i++) {
        r.coeffRef(i, 0) = integral(betta, f, x_from, x_to, i, 0, h);
        for(int j = 0; j <= n; j++) {
            M.coeffRef(i, j) = -lambda*integral(betta, alpha, x_from, x_to, i, j, h);
        }
    }
    for(int i = 0; i <= n; i++)
        M.coeffRef(i, i) = 1 + M.coeffRef(i, i);

    Matrix c = M.inverse()*r;

    y.clear();
    for(int k = 0; k < abs(x_to-x_from)/h; k++) {
        double x = x_from + sign(x_to-x_from)*k*h;
        double res = 0;
        for (int i = 0; i <= n; i++)
            res += c.coeffRef(i, 0)*alpha(x, i);
        res = res*lambda + f(x, 0);

        y.push(x, res);
    }
}

void test_Degenerate_Fred_2(ComplexPlot& y, double h = 0.01) {
    y.clear();
    for(int k = 0; k < 0.5/h; k++) {
        double x = k*h;
        y.push(x, 1.993 - 0.0833*x*x + 0.0007*x*x*x*x);
    }
}

int main() {

    double x_from = 0;
    double x_to = 0.5;
    double lambda = 1;

    double h = 0.001;

    ComplexPlot y1;
    ComplexPlot y2;
    ComplexPlot y3;
    ComplexPlot y4;
    ComplexPlot y5;

    ComplexPlot y_analytical;

    Degenerate_Fred_2(y1, x_from, x_to, lambda, &alpha_i, &betta_i, &funct, 1, h);
    Degenerate_Fred_2(y2, x_from, x_to, lambda, &alpha_i, &betta_i, &funct, 2, h);
    Degenerate_Fred_2(y3, x_from, x_to, lambda, &alpha_i, &betta_i, &funct, 3, h);
    Degenerate_Fred_2(y4, x_from, x_to, lambda, &alpha_i, &betta_i, &funct, 4, h);
    Degenerate_Fred_2(y5, x_from, x_to, lambda, &alpha_i, &betta_i, &funct, 5, h);

    test_Degenerate_Fred_2(y_analytical, h);

    saveVectorPoint2DToFile(y1.real(), "y1.dat");
    saveVectorPoint2DToFile(y2.real(), "y2.dat");
    saveVectorPoint2DToFile(y3.real(), "y3.dat");
    saveVectorPoint2DToFile(y4.real(), "y4.dat");
    saveVectorPoint2DToFile(y5.real(), "y5.dat");

    saveVectorPoint2DToFile(y_analytical.real(), "y_analytical.dat");

    GnuplotPipe gp;

    gp.sendLine(R"(plot "y1.dat" with lines, "y_analytical.dat" with lines)");

    //gp.sendLine(R"(plot "y1.dat" with lines, "y2.dat" with lines, "y3.dat" with lines, "y4.dat" with lines, "y5.dat" with lines,)");
}