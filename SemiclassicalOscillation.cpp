//
// Created by ivan- on 31.03.2020.
//
#include <cmath>
#include <iostream>
#include <functional>

using namespace std;

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

double integral(double (*pFunction)(double x), double from, double to, double h = 0.01) {
    double d_summ = 0;
    double d_x = from;
    while(d_x < to) {
        d_summ += pFunction(d_x) + 4*pFunction(d_x + h) + pFunction(d_x + 2*h);
        d_x += 2*h;
    }
    return d_summ * h / 3;
}



double En (int n, double E, double gamma) {
    double d_from = x_in(E);
    double d_to   = x_out(E);

    double d_En = E_ERROR;
    double d_F = 0;
    do {
        d_F = gamma * integral( [](double x) -> double {
            return 1 - V(x);
            }, d_from, d_to);

    } while (abs(d_F) > E_ERROR);

}

int main() {
    //cout << integral(&testFunc, 0, 10) << endl;


}