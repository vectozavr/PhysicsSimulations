//
// Created by Иван Ильин on 17.09.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;

double U_0 = 1.0f;
double H = 1.0f;
double m = 1.0f;
double A = 1.0f;

double ctg(double x) {
    return cos(x) / sin(x);
}

double f(double ksi) {
    return ctg(sqrt(2*(1-ksi)*m*A*A*U_0/(H*H))) - sqrt(-1 + 1/ksi);
}

double df(double ksi, double h, double fun(double ksi)) {
    return (fun(ksi+h) - fun(ksi-h))/(2*h);
}

int main() {
    double eps = 0.0000001f; // 10^-8

    double a = 0.1f, b = 0.9f; //study area

    double a1 = 0.1f, b1 = 0.0f;
    int n = 0;

    // Дихотомия
    while (abs(b - a) > eps) {
        n++;
        a1 = (a + b) / 2.0f;
        b1 = (a + b) / 2.0f;
        if (f(a1) * f(b) <= 0)
            a = a1;
        else
            b = b1;
    }

    // Метод простых итераций
    double ksi = 0, ksi1 = 0.2, la = 0.1;
    int k = 0;

    while ( abs(ksi - ksi1) > eps ){
        ksi = ksi1;
        ksi1 = ksi - la * f(ksi);
        k++;
    }

    // Метод Ньютона
    double ksiN = 0.0f, ksi1N = 0.5064675f, h = 0.0001f;
    int s = 0;

    while(abs(ksiN - ksi1N) > eps) {
        ksiN = ksi1N;
        ksi1N = ksiN - f(ksiN) / df(ksiN, h, f);
        s++;
    }

    cout << "Dichotomy method   (n = "    << n << ")    E = " << a       << endl;
    cout << "Simple iteration   (n = "    << k << ")    E = " << ksi     << endl;
    cout << "Newton method      (n = "       << s << ")     E = " << ksiN    << endl;
}